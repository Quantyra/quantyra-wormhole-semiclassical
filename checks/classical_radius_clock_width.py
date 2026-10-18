"""S140 classical-radius line search and clock-width test.

The checker calls one calculation on the retained horizon, disk,
same-source, and seam notes. It does not assign b0, and it does not
evaluate an area budget.
"""
from __future__ import annotations
import base64, contextlib, hashlib, json, platform, re, sys, traceback
from datetime import datetime, timezone
from pathlib import Path
import sympy as S

ROOT = Path(__file__).resolve().parents[1]
RECEIPT = ROOT / "evidence/classical-radius-clock-width.json"
LOG = ROOT / "evidence/classical-radius-clock-width.log"
INPUTS = [
    "checks/classical_radius_clock_width.py",
    "notes/classical-radius-clock-width.md",
    "notes/stationary-horizon-product-matching.md",
    "notes/stationary-disk-and-normal-inverse.md",
    "notes/stationary-same-source-boundary-charge.md",
    "notes/stationary-variable-radius-seam.md",
]
CONTRACT = {
    "acceptance": "the identity, the strict lower bound, the horizon line search, and the clock width agree",
    "numerical_tolerances": "none",
    "scope": "one frozen spectral point and the retained horizon equations",
    "ordinary_arguments": [
        "comparison, hardest step, observable, and test recorded before the fit",
        "line search for one classical extremal radius",
        "width of kappa_M from the recorded bounds",
    ],
    "unresolved": [
        "one numerical classical extremal radius b0",
        "S074.2 freeze",
        "area budget",
    ],
    "not_set_to_zero": ["kappa", "lambda", "G", "m_num", "h_R"],
    "area_budget": "not evaluated; no replacement target was written before the test",
    "failure_policy": "append initialized scientific runs and retain failed raw input bytes",
    "hash_policy": "SHA256 of raw UTF8/LF/noBOM bytes without normalization",
}
ROWS = []
R = S.Rational
ASSIGNMENT = re.compile(
    r"(?<![A-Za-z0-9_])(b0|b_classical)\s*=\s*"
    r"([0-9]+(?:/[0-9]+)?(?:\.[0-9]+)?(?:[eE][+-]?[0-9]+)?)"
)
ZERO_ASSIGNMENT = {
    "kappa": re.compile(r"(?<![A-Za-z0-9_])kappa\s*=\s*0\b"),
    "m_num": re.compile(r"(?<![A-Za-z0-9_])m_num\s*=\s*0\b"),
    "h_R": re.compile(r"\|\|h_R\|\|\s*=\s*0\b"),
    "lambda": re.compile(r"(?<![A-Za-z0-9_])lambda\s*=\s*0\b"),
    "G": re.compile(r"(?<![A-Za-z0-9_])G\s*=\s*0\b"),
}


def check(label, condition, detail=""):
    if any(row["label"] == label for row in ROWS):
        raise ValueError("duplicate check label: " + label)
    row = {"label": label, "passed": bool(condition), "detail": str(detail)}
    ROWS.append(row)
    print(json.dumps(row, sort_keys=True), flush=True)


def zero(label, expression):
    expression = S.sympify(expression)
    if expression.atoms(S.Float):
        raise ValueError("inexact Float atom in " + label)
    value = S.simplify(S.together(expression))
    check(label, value == 0, value)


def line_search(text):
    """Return numeral assignments of b0 or b_classical. Inequalities are excluded."""
    found = []
    for lineno, line in enumerate(text.splitlines(), 1):
        if "<" in line or ">" in line:
            continue
        if "mu0=1/b0" in line or "mu0=1/b_classical" in line:
            continue
        for match in ASSIGNMENT.finditer(line):
            found.append({"line": lineno, "name": match.group(1), "value": match.group(2)})
    return found


def evaluate(horizon, disk, same, seam):
    """Pure S140 calculation. The texts are the retained notes."""
    if horizon.find("mu0=1/b0") < 0:
        raise ValueError("horizon scale missing")
    if horizon.find("not reset to the corrected radius") < 0:
        raise ValueError("corrected-radius sentence missing")
    if seam.find("G=1,  lambda=(1999/2000)*2000^{-10}") < 0:
        raise ValueError("frozen point missing")
    disk_bound = re.search(r"k=(\d+),\s*lambda b_classical\^2>k\^3/4", disk)
    if disk_bound is None:
        raise ValueError("disk lower bound missing")
    if same.find("1<lambda b_classical^2<2k^3") < 0:
        raise ValueError("same-source interval missing")
    k = S.Integer(int(disk_bound.group(1)))
    G = S.Integer(1)
    lam = R(1999, 2000**11)
    if lam != R(1999, 2000) * S.Integer(2000) ** (-10):
        raise ValueError("frozen lambda does not match 1999/2000^{11}")
    ell_squared = S.Integer(1) / lam
    lower_product = k**3 / 4
    b0_sq_endpoint = lower_product / lam
    assignments = line_search(horizon) + line_search(disk) + line_search(same)
    endpoint_used = any(item["value"] == str(b0_sq_endpoint) for item in assignments)
    width_zero = len(assignments) == 1 and not endpoint_used
    disk_width = R(54, 5) - R(39, 5)
    bound_ratio = (2 * k**3) / lower_product
    log_width = S.log(bound_ratio) / 2
    if any(expr.atoms(S.Float) for expr in (lam, ell_squared, b0_sq_endpoint, log_width)):
        raise ValueError("inexact Float atom in the radius calculation")
    return {
        "G": G,
        "lambda": lam,
        "ell_squared": ell_squared,
        "k": k,
        "lower_product": lower_product,
        "b0_sq_endpoint": b0_sq_endpoint,
        "assignments": assignments,
        "endpoint_used": endpoint_used,
        "b0": assignments[0]["value"] if width_zero else None,
        "width_zero": width_zero,
        "width_positive": not width_zero,
        "disk_enclosure_width": disk_width,
        "bound_ratio": bound_ratio,
        "log_width": log_width,
        "missing": None if width_zero else "numerical b0",
        "classification": "one numerical b0" if width_zero else "unavailable physical input",
        "kappa_is_zero": False,
        "m_num_is_zero": False,
        "h_R_norm_is_zero": False,
        "area_budget": "not evaluated",
    }


def classical_radius_clock_width(horizon, disk, same, seam, argument):
    result = evaluate(horizon, disk, same, seam)
    b0, lam_s = S.symbols("b0 lam", positive=True)
    mu0 = 1 / b0
    ell = lam_s ** (-R(1, 2))
    zero("comparison_identity", S.log(mu0 * ell) + S.log(b0 * S.sqrt(lam_s)))
    zero("clock_shift_identity", (-R(6, 5) - S.log(mu0 * ell)) - (-R(6, 5) + S.log(b0 * S.sqrt(lam_s))))
    zero("ell_squared_times_lambda", result["ell_squared"] * result["lambda"] - 1)
    zero("frozen_lambda", result["lambda"] - R(1999, 2000**11))
    check("frozen_G_is_one", result["G"] == 1 and result["G"] != 0, result["G"])
    check("lambda_not_zero", result["lambda"] != 0, result["lambda"])
    zero("lower_product_endpoint", result["lambda"] * result["b0_sq_endpoint"] - result["k"] ** 3 / 4)
    check("lower_product_is_2000_cubed_over_4", result["lower_product"] == R(2000**3, 4), result["lower_product"])
    check("lower_bound_is_strict", "lambda b_classical^2>k^3/4" in disk and "lambda b_classical^2>=" not in disk)
    check("integer_guard_below_lower_product", 3**18 < result["lower_product"], 3**18)
    check("endpoint_not_used_as_b0", result["endpoint_used"] is False and result["b0"] is None)
    check("corrected_radius_not_substituted", "b0=b_star" not in horizon.replace(" ", "") and "not reset to the corrected radius" in horizon)
    check("horizon_equations_present", "8piG pA,h" in horizon and "kappa=exp(alpha_h)" in horizon)
    check("no_numerical_b0", result["assignments"] == [], result["assignments"])
    check("width_zero_only_with_one_b0", result["width_zero"] == (len(result["assignments"]) == 1 and not result["endpoint_used"]))
    check("width_is_not_zero", result["width_zero"] is False)
    check("width_positive", result["width_positive"] is True)
    zero("disk_enclosure_width", result["disk_enclosure_width"] - 3)
    zero("disk_lower_shift", R(39, 5) - (9 - R(6, 5)))
    zero("disk_upper_shift", R(54, 5) - (12 - R(6, 5)))
    check("disk_enclosure_width_positive", result["disk_enclosure_width"] == 3 and result["disk_enclosure_width"] != 0)
    zero("combined_bound_ratio", result["bound_ratio"] - 8)
    check("combined_bounds_do_not_meet", result["bound_ratio"] > 1, result["bound_ratio"])
    zero("combined_log_width", result["log_width"] - R(3, 2) * S.log(2))
    check("combined_log_width_positive", result["log_width"] > 0, result["log_width"])
    check("same_source_upper_bound_recorded", "1<lambda b_classical^2<2k^3" in same)
    check("log_bound_integer_guard", 8**24 > 2 * int(result["k"]) ** 3 * 3**24)
    check("missing_input_is_numerical_b0", result["missing"] == "numerical b0" and result["classification"] == "unavailable physical input")
    check("kappa_not_set_to_zero", result["kappa_is_zero"] is False and ZERO_ASSIGNMENT["kappa"].search(argument) is None)
    check("m_num_not_set_to_zero", result["m_num_is_zero"] is False and ZERO_ASSIGNMENT["m_num"].search(argument) is None)
    check("h_R_norm_not_set_to_zero", result["h_R_norm_is_zero"] is False and ZERO_ASSIGNMENT["h_R"].search(argument) is None)
    check("lambda_assignment_not_zero", ZERO_ASSIGNMENT["lambda"].search(argument) is None)
    check("G_assignment_not_zero", ZERO_ASSIGNMENT["G"].search(argument) is None)
    protocol = argument.find("## 1. Comparison, hardest step, observable, and test")
    verdict = argument.find("## 5. Classification")
    check("protocol_before_verdict", 0 <= protocol < verdict)
    head = argument[protocol:verdict]
    check("protocol_states_comparison", "log(mu0 ell)=-log(b0 sqrt(lambda))" in head)
    check("protocol_states_hardest_step", "classical extremal radius" in head and "2000^3/4" in head)
    check("protocol_states_observable", "width of" in head and "kappa_M=-6/5-log(mu0 ell)" in head)
    check("protocol_states_test", "line search" in head)
    tail = argument[verdict:]
    check("verdict_matches_search", "unavailable physical input" in tail and "width is not 0" in tail and "numerical classical extremal radius" in tail)
    check("area_budget_unevaluated", result["area_budget"] == "not evaluated" and "area budget is unevaluated" in tail and "No replacement area target" in tail)
    return result


def hashes(data):
    return {name: hashlib.sha256(raw).hexdigest() for name, raw in data.items()}


def main():
    before = {name: (ROOT / name).read_bytes() for name in INPUTS}
    existing = json.loads(RECEIPT.read_text(encoding="utf-8")) if RECEIPT.exists() else {"runs": []}
    run = {
        "started_utc": datetime.now(timezone.utc).isoformat(),
        "contract": CONTRACT,
        "command": sys.orig_argv,
        "dont_write_bytecode": sys.dont_write_bytecode,
        "cwd": str(Path.cwd()),
        "python": platform.python_version(),
        "sympy": S.__version__,
        "execution_inputs_sha256": hashes(before),
    }
    with LOG.open("a", encoding="utf-8", newline="\n") as log, contextlib.redirect_stdout(log), contextlib.redirect_stderr(log):
        print("RUN " + json.dumps(run, sort_keys=True), flush=True)
        try:
            check("actual_bytecode_disabled", sys.dont_write_bytecode and "-B" in sys.orig_argv)
            check(
                "inputs_strict_UTF8_LF",
                all(
                    raw.decode("utf-8") is not None and b"\r" not in raw and not raw.startswith(b"\xef\xbb\xbf")
                    for raw in before.values()
                ),
            )
            texts = {name: raw.decode("utf-8") for name, raw in before.items()}
            classical_radius_clock_width(
                texts["notes/stationary-horizon-product-matching.md"],
                texts["notes/stationary-disk-and-normal-inverse.md"],
                texts["notes/stationary-same-source-boundary-charge.md"],
                texts["notes/stationary-variable-radius-seam.md"],
                texts["notes/classical-radius-clock-width.md"],
            )
        except Exception:
            traceback.print_exc()
            check("uncaught_exception", False)
        run["final_inputs_sha256"] = hashes({name: (ROOT / name).read_bytes() for name in INPUTS})
        print("FINAL_INPUTS " + json.dumps(run["final_inputs_sha256"], sort_keys=True), flush=True)
        check("frozen_inputs_unchanged", run["execution_inputs_sha256"] == run["final_inputs_sha256"])
        run["checks"] = ROWS
        run["failed_checks"] = [row for row in ROWS if not row["passed"]]
        code = int(bool(run["failed_checks"]))
        run["status"] = "failed" if code else "passed"
        run["terminal_exit_status"] = code
        run["finished_utc"] = datetime.now(timezone.utc).isoformat()
        if code:
            run["failed_inputs_base64"] = {name: base64.b64encode(raw).decode("ascii") for name, raw in before.items()}
        print(
            "TERMINAL "
            + json.dumps(
                {
                    "status": run["status"],
                    "checks": len(ROWS),
                    "failures": len(run["failed_checks"]),
                    "terminal_exit_status": code,
                    "finished_utc": run["finished_utc"],
                },
                sort_keys=True,
            ),
            flush=True,
        )
    existing["runs"].append(run)
    RECEIPT.write_text(json.dumps(existing, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    print(json.dumps({"status": run["status"], "checks": len(ROWS), "failures": len(run["failed_checks"]), "receipt": str(RECEIPT)}, sort_keys=True))
    return code


if __name__ == "__main__":
    sys.exit(main())
