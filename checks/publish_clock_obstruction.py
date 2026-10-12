"""Publish the clock-obstruction preprint as its own Zenodo record.

--compile runs the portable MiKTeX route on the manuscript.
--reserve is the only new-record creation route.
--prepare builds the local source archive after the reserved DOI is in the TeX.
--publish --commit HASH uploads and publishes only that reserved record.
--verify-public is read-only.
No Lean or Lake calls. Credentials stay in memory. No automatic POST retry.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import zipfile
from urllib.parse import urlparse

import requests

ROOT = Path(__file__).resolve().parents[1]
NULL_CONE_CHECKS = ROOT.parent / "quantyra-null-cone-lean" / "checks"
PAPER = ROOT / "manuscript" / "clock-mode-obstruction"
PACKAGE = PAPER / "deposit" / "v0.1.0"
TEX = PAPER / "clock-mode-obstruction.tex"
PDF_NAME = "clock-mode-obstruction.pdf"
STEM = "clock-mode-obstruction"
API = "https://zenodo.org/api"
TITLE = "The symmetric MMP mouth clock is not an independent boundary mode"
PROOF = "e02e76a99b06a20c305d71dcc96117530cce941b"
OLD_FAMILIES = {"23206772", "23249815", "23265066", "23202762"}
PDFTOTOTEXT = Path(
    r"C:\Users\dfred\.quantyra\tools\miktex-draft-20261009"
    r"\portable\texmfs\install\miktex\bin\x64\pdftotext.exe"
)
PDFLATEX = Path(
    r"C:\Users\dfred\.quantyra\tools\miktex-draft-20261009"
    r"\portable\texmfs\install\miktex\bin\x64\pdflatex.exe"
)
BOUNDARY = (
    "already known",
    "No shortcut, formation, propulsion, occupancy, TRL, launch, seed, or revenue result is claimed",
)


def sha(data):
    return hashlib.sha256(data).hexdigest()


def save(path, value):
    path.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8", newline="\n")


def read(path):
    return json.loads(path.read_text(encoding="utf-8"))


def trusted(url):
    parsed = urlparse(url)
    assert parsed.scheme == "https" and parsed.netloc == "zenodo.org"
    assert parsed.path.startswith("/api/") and not parsed.query and not parsed.fragment
    return url


def validate_metadata(metadata):
    assert metadata["title"] == TITLE
    assert metadata["version"] == "0.1.0"
    assert metadata["upload_type"] == "publication" and metadata["publication_type"] == "preprint"
    assert metadata["license"] == "cc-by-4.0" and metadata["access_right"] == "open"
    assert metadata["publication_date"] == "2026-10-12"
    assert metadata["creators"] == [{"name": "Fredriksen, Daniel Eric", "affiliation": "Quantyra Inc"}]
    assert "doi" not in metadata
    description = metadata["description"]
    for phrase in BOUNDARY:
        assert phrase in description
    assert "new theorem" in description
    supplement = (
        "https://github.com/Quantyra/quantyra-wormhole-semiclassical/tree/" + PROOF
    )
    assert any(
        item["identifier"] == supplement and item["relation"] == "isSupplementedBy"
        for item in metadata["related_identifiers"]
    )


def validate_identity(record, reservation=None):
    assert str(record["conceptrecid"]) not in OLD_FAMILIES
    assert record["metadata"]["title"] == TITLE
    if reservation:
        assert record["id"] == reservation["id"]
        assert str(record["conceptrecid"]) == reservation["conceptrecid"]
        doi = record["metadata"].get("doi") or record["metadata"]["prereserve_doi"]["doi"]
        assert doi == reservation["doi"]


def source_text():
    return TEX.read_text(encoding="utf-8")


def pdf_text(pdf):
    if not PDFTOTOTEXT.is_file():
        raise SystemExit(f"pdftotext is not at the satellite TeX route: {PDFTOTOTEXT}")
    result = subprocess.run(
        [str(PDFTOTOTEXT), "-layout", str(pdf), "-"],
        capture_output=True,
        check=False,
    )
    if result.returncode != 0:
        raise SystemExit("pdftotext failed")
    return result.stdout.decode("utf-8", errors="replace")


def assert_claim_boundary(text):
    collapsed = " ".join(text.split())
    assert "The symmetric MMP mouth clock is not an independent boundary mode" in collapsed
    for phrase in BOUNDARY:
        assert phrase in collapsed
    assert "Publication is not selected" not in collapsed


def compile_pdf():
    if not PDFLATEX.is_file():
        raise SystemExit(f"pdflatex is not at the satellite TeX route: {PDFLATEX}")
    PACKAGE.mkdir(parents=True, exist_ok=True)
    log_path = PACKAGE / "pdflatex.log"
    for _pass in (1, 2):
        result = subprocess.run(
            [
                str(PDFLATEX),
                "-interaction=nonstopmode",
                "-halt-on-error",
                f"-output-directory={PACKAGE}",
                TEX.name,
            ],
            cwd=PAPER,
            capture_output=True,
        )
        log_path.write_bytes(result.stdout + b"\n" + result.stderr)
        if result.returncode != 0:
            raise SystemExit(f"pdflatex failed on pass {_pass}; log retained at {log_path}")
    pdf = PACKAGE / PDF_NAME
    if not pdf.is_file() or not pdf.read_bytes().startswith(b"%PDF-"):
        raise SystemExit("pdflatex did not write a PDF")
    log = log_path.read_text(encoding="utf-8", errors="replace")
    if "Output written on" not in log:
        raise SystemExit("pdflatex log has no output line")
    for name in ("aux", "out"):
        extra = PACKAGE / f"{STEM}.{name}"
        if extra.exists():
            extra.unlink()
    print(json.dumps({"pdf": str(pdf), "bytes": pdf.stat().st_size, "sha256": sha(pdf.read_bytes())}))


def authenticated_request():
    sys.path.insert(0, str(NULL_CONE_CHECKS))
    from zenodo_credentials import aws_token

    token, _metadata = aws_token()
    session = requests.Session()
    session.headers["Authorization"] = "Bearer " + token

    def request(method, url, **kwargs):
        try:
            response = session.request(
                method, trusted(url), timeout=90, allow_redirects=False, **kwargs
            )
        except requests.RequestException:
            raise RuntimeError("Zenodo transport failure; outcome may be unknown; no automatic retry") from None
        if not 200 <= response.status_code < 300:
            raise RuntimeError(f"Zenodo {method} HTTP {response.status_code}; body withheld")
        return response

    return request


def reservation_path():
    return PACKAGE / "draft-record.json"


def intent_path():
    return PACKAGE / "creation-intent.json"


def reserve(metadata):
    request = authenticated_request()
    path = reservation_path()
    if path.exists():
        reservation = read(path)
        draft = request("GET", f"{API}/deposit/depositions/{reservation['id']}").json()
        validate_identity(draft, reservation)
        print(json.dumps(reservation, indent=2))
        return reservation
    candidates = []
    for page in range(1, 51):
        deposits = request("GET", f"{API}/deposit/depositions", params={"size": 100, "page": page}).json()
        assert isinstance(deposits, list)
        candidates.extend(item for item in deposits if item.get("metadata", {}).get("title") == TITLE)
        if len(deposits) < 100:
            break
    else:
        raise RuntimeError("Deposit listing incomplete; no creation attempted")
    assert len(candidates) <= 1, "Multiple matching records require inspection"
    if candidates:
        draft = request("GET", f"{API}/deposit/depositions/{candidates[0]['id']}").json()
    else:
        assert not intent_path().exists(), "Unresolved creation intent; inspect the account before retrying"
        PACKAGE.mkdir(parents=True, exist_ok=True)
        save(intent_path(), {"title": TITLE, "created_utc": datetime.now(timezone.utc).isoformat()})
        draft = request("POST", f"{API}/deposit/depositions", json={"metadata": metadata}).json()
    validate_identity(draft)
    assert not draft["submitted"], "Existing publication requires inspection, not a duplicate"
    reservation = {
        "id": draft["id"],
        "doi": draft["metadata"]["prereserve_doi"]["doi"],
        "conceptrecid": str(draft["conceptrecid"]),
        "conceptdoi": "10.5281/zenodo." + str(draft["conceptrecid"]),
        "independent_manuscript_family": True,
        "version": "0.1.0",
    }
    save(path, reservation)
    print(json.dumps(reservation, indent=2))
    return reservation


def stamp_doi(doi):
    assert doi.startswith("10.5281/zenodo.") and doi.removeprefix("10.5281/zenodo.").isdigit()
    text = source_text()
    token = "10.5281/zenodo.PENDING"
    if doi in text and token not in text:
        return
    assert text.count(token) == 1
    TEX.write_text(text.replace(token, doi, 1), encoding="utf-8", newline="\n")


def license_file():
    destination = PACKAGE / "LICENSE-CC-BY-4.0.txt"
    source = NULL_CONE_CHECKS.parent / "LICENSES" / "CC-BY-4.0.txt"
    if source.is_file():
        shutil.copyfile(source, destination)
    elif not destination.is_file():
        destination.write_text(
            "This manuscript and source bundle are licensed under CC-BY-4.0.\n"
            "https://creativecommons.org/licenses/by/4.0/\n",
            encoding="utf-8",
            newline="\n",
        )
    return destination


def prepare(metadata):
    assert not (PACKAGE / "published-record.json").exists(), "Published package is immutable"
    reservation = read(reservation_path())
    stamp_doi(reservation["doi"])
    text = source_text()
    assert_claim_boundary(text)
    assert reservation["doi"] in text
    pdf = PACKAGE / PDF_NAME
    rendered = pdf_text(pdf)
    assert pdf.is_file() and reservation["doi"] in rendered
    assert_claim_boundary(rendered)
    license_path = license_file()
    payload = [TEX, PAPER / "claim-to-evidence.md", pdf, PACKAGE / "zenodo-metadata.json", license_path]
    manifest = {
        "title": TITLE,
        "version": "0.1.0",
        "doi": reservation["doi"],
        "proof_commit": PROOF,
        "metadata_sha256": sha((PACKAGE / "zenodo-metadata.json").read_bytes()),
        "files": {
            path.relative_to(ROOT).as_posix(): {"sha256": sha(path.read_bytes()), "bytes": path.stat().st_size}
            for path in payload
        },
    }
    save(PACKAGE / "manifest.json", manifest)
    archive_path = PACKAGE / f"{STEM}-v0.1.0-source.zip"
    with zipfile.ZipFile(archive_path, "w") as archive:
        for path in payload + [PACKAGE / "manifest.json"]:
            entry = zipfile.ZipInfo(path.relative_to(ROOT).as_posix(), (2026, 10, 12, 0, 0, 0))
            entry.compress_type = zipfile.ZIP_DEFLATED
            entry.external_attr = 0o100644 << 16
            archive.writestr(entry, path.read_bytes())
    with zipfile.ZipFile(archive_path) as archive:
        assert archive.testzip() is None
        for name, info in manifest["files"].items():
            assert sha(archive.read(name)) == info["sha256"]
    save(
        PACKAGE / "bundle.json",
        {
            "file": archive_path.name,
            "sha256": sha(archive_path.read_bytes()),
            "bytes": archive_path.stat().st_size,
            "members": len(payload) + 1,
        },
    )
    print("PASS: claim boundary, reserved DOI and deterministic source archive")


def frozen_bytes(commit, path):
    return subprocess.check_output(["git", "show", f"{commit}:{path.relative_to(ROOT).as_posix()}"], cwd=ROOT)


def frozen_payload(commit):
    assert len(commit) == 40 and all(c in "0123456789abcdef" for c in commit)
    for name in ("manifest.json", "bundle.json", "draft-record.json", "zenodo-metadata.json"):
        assert frozen_bytes(commit, PACKAGE / name) == (PACKAGE / name).read_bytes(), name
    metadata = read(PACKAGE / "zenodo-metadata.json")
    validate_metadata(metadata)
    manifest = read(PACKAGE / "manifest.json")
    assert manifest["proof_commit"] == PROOF
    assert manifest["doi"] == read(reservation_path())["doi"]
    for name, info in manifest["files"].items():
        path = ROOT / name
        assert sha(frozen_bytes(commit, path)) == info["sha256"] == sha(path.read_bytes()), name
    bundle = read(PACKAGE / "bundle.json")
    archive = PACKAGE / bundle["file"]
    assert sha(frozen_bytes(commit, archive)) == bundle["sha256"] == sha(archive.read_bytes())
    pdf = PACKAGE / PDF_NAME
    assert sha(frozen_bytes(commit, pdf)) == sha(pdf.read_bytes())
    print("PASS: frozen Git metadata, manifest and upload bytes", flush=True)
    return metadata, read(reservation_path()), [pdf, archive]


def public_record(record_id):
    response = requests.get(f"{API}/records/{record_id}", timeout=60)
    response.raise_for_status()
    return response.json()


def public_files(record):
    result = {}
    for entry in record["files"]:
        response = requests.get(trusted(entry["links"]["self"]), timeout=90)
        response.raise_for_status()
        assert entry["checksum"] == "md5:" + hashlib.md5(response.content).hexdigest()
        result[entry["key"]] = {"sha256": sha(response.content), "zenodo_checksum": entry["checksum"]}
    return result


def verify_public(commit, metadata, reservation, payload):
    record = public_record(reservation["id"])
    validate_identity(record, reservation)
    actual = record["metadata"]
    for key in ("title", "version", "creators", "publication_date", "description", "access_right", "keywords"):
        assert actual[key] == metadata[key], key
    assert actual["license"]["id"].lower() == "cc-by-4.0"
    assert actual["resource_type"]["type"] == "publication"
    assert actual["resource_type"]["subtype"] == "preprint"
    for relation in metadata["related_identifiers"]:
        assert relation in actual["related_identifiers"], relation
    verified = public_files(record)
    assert set(verified) == {path.name for path in payload}
    for path in payload:
        assert verified[path.name]["sha256"] == sha(path.read_bytes())
    doi = actual["doi"]
    resolution = requests.get("https://doi.org/" + doi, timeout=60, allow_redirects=True)
    assert resolution.status_code == 200
    pdf_url = f"https://zenodo.org/records/{record['id']}/files/{PDF_NAME}?download=1"
    downloaded = requests.get(pdf_url, timeout=90)
    assert downloaded.status_code == 200
    assert sha(downloaded.content) == sha(payload[0].read_bytes())
    zip_url = f"https://zenodo.org/records/{record['id']}/files/{payload[1].name}?download=1"
    downloaded_zip = requests.get(zip_url, timeout=90)
    assert downloaded_zip.status_code == 200
    assert sha(downloaded_zip.content) == sha(payload[1].read_bytes())
    receipt = {
        "id": record["id"],
        "doi": doi,
        "conceptdoi": record["conceptdoi"],
        "url": f"https://zenodo.org/records/{record['id']}",
        "version": "0.1.0",
        "source_commit": commit,
        "proof_commit": PROOF,
        "independent_manuscript_family_verified": True,
        "verified_utc": datetime.now(timezone.utc).isoformat(),
        "files": verified,
        "pdf_download_sha256": sha(downloaded.content),
        "source_download_sha256": sha(downloaded_zip.content),
    }
    save(PACKAGE / "public-record.json", record)
    save(PACKAGE / "published-record.json", receipt)
    save(
        PACKAGE / "resolution-check.json",
        {
            "doi": doi,
            "resolved_url": resolution.url,
            "doi_http_status": resolution.status_code,
            "pdf_url": pdf_url,
            "pdf_http_status": downloaded.status_code,
            "pdf_sha256": sha(downloaded.content),
            "source_url": zip_url,
            "source_http_status": downloaded_zip.status_code,
            "source_sha256": sha(downloaded_zip.content),
            "verified_utc": receipt["verified_utc"],
        },
    )
    print(json.dumps({"doi": doi, "url": receipt["url"], "files": verified}, indent=2))


def publish(commit, metadata, reservation, payload):
    request = authenticated_request()
    record_url = f"{API}/deposit/depositions/{reservation['id']}"
    draft = request("GET", record_url).json()
    validate_identity(draft, reservation)
    if not draft["submitted"]:
        existing = {entry["filename"]: entry for entry in draft["files"]}
        assert set(existing) <= {path.name for path in payload}, "Unexpected files: no changes made"
        for path in payload:
            if path.name in existing:
                entry = existing[path.name]
                if entry["checksum"].removeprefix("md5:") == hashlib.md5(path.read_bytes()).hexdigest():
                    continue
                request("DELETE", record_url + "/files/" + str(entry["id"]))
            with path.open("rb") as stream:
                request("POST", record_url + "/files", data={"name": path.name}, files={"file": (path.name, stream)})
        draft = request("PUT", record_url, json={"metadata": metadata}).json()
        validate_identity(draft, reservation)
        files = {entry["filename"]: entry for entry in draft["files"]}
        assert set(files) == {path.name for path in payload}
        for path in payload:
            assert files[path.name]["checksum"].removeprefix("md5:") == hashlib.md5(path.read_bytes()).hexdigest()
        for key in ("title", "version", "creators", "publication_date", "description", "access_right", "keywords"):
            assert draft["metadata"][key] == metadata[key], key
        assert draft["metadata"]["license"].lower() == "cc-by-4.0"
        print("PASS: reserved draft, uploaded checksums and metadata; publishing", flush=True)
        request("POST", record_url + "/actions/publish")
    verify_public(commit, metadata, reservation, payload)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group(required=True)
    for mode in ("compile", "reserve", "prepare", "publish", "verify-public"):
        modes.add_argument("--" + mode, action="store_true")
    parser.add_argument("--commit")
    args = parser.parse_args(argv)
    metadata = read(PACKAGE / "zenodo-metadata.json")
    validate_metadata(metadata)
    if args.compile:
        compile_pdf()
        return 0
    if args.reserve:
        assert not args.commit
        reserve(metadata)
        return 0
    if args.prepare:
        assert not args.commit
        prepare(metadata)
        return 0
    assert args.commit, "--commit is required"
    metadata, reservation, payload = frozen_payload(args.commit)
    if args.verify_public:
        verify_public(args.commit, metadata, reservation, payload)
        return 0
    publish(args.commit, metadata, reservation, payload)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
