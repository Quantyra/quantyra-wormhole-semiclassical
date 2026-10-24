"""Read research PDFs in memory and retain source bytes/text only in private S3.

Requires boto3 and pypdf. No PDF or extracted full text is written locally.
"""
import argparse
import hashlib
import io
import json
from pathlib import Path
import subprocess
import sys
import zipfile

import boto3
from pypdf import PdfReader

BUCKET = 'quantyra-research-archive-063280428495-us-east-1'
OWNER = '063280428495'
PREFIX = 'uap-release-analysis/d153-assessment/2026-10-24'
ROOT = Path(__file__).resolve().parents[1]
INDEX = ROOT / 'evidence/d153-source-custody.json'


class S3Range(io.RawIOBase):
    def __init__(self, s3, key, version, size):
        self.s3, self.key, self.version, self.size, self.pos = s3, key, version, size, 0

    def seekable(self):
        return True

    def tell(self):
        return self.pos

    def seek(self, offset, whence=0):
        self.pos = offset if whence == 0 else (self.pos if whence == 1 else self.size) + offset
        return self.pos

    def read(self, size=-1):
        end = self.size if size < 0 else min(self.size, self.pos + size)
        if end <= self.pos:
            return b''
        response = self.s3.get_object(Bucket=BUCKET, Key=self.key, VersionId=self.version,
            ExpectedBucketOwner=OWNER, Range=f'bytes={self.pos}-{end-1}')
        data = response['Body'].read()
        self.pos += len(data)
        return data


def archive(s3, name, pdf, provenance):
    assert pdf.startswith(b'%PDF-'), 'Not a PDF'
    pages = [p.extract_text() or '' for p in PdfReader(io.BytesIO(pdf)).pages]
    extracted = json.dumps(pages, ensure_ascii=False).encode('utf-8')
    saved = {}
    for suffix, data, content_type in [('pdf', pdf, 'application/pdf'), ('pages.json', extracted, 'application/json')]:
        key = f'{PREFIX}/{name}.{suffix}'
        result = s3.put_object(Bucket=BUCKET, Key=key, Body=data, ContentType=content_type,
            ServerSideEncryption='AES256', ExpectedBucketOwner=OWNER)
        version = result['VersionId']
        check = s3.get_object(Bucket=BUCKET, Key=key, VersionId=version,
            ExpectedBucketOwner=OWNER)['Body'].read()
        digest = hashlib.sha256(data).hexdigest()
        assert hashlib.sha256(check).hexdigest() == digest
        saved[suffix] = {'key': key, 'version_id': version, 'bytes': len(data), 'sha256': digest}
    index = json.loads(INDEX.read_text()) if INDEX.exists() else {'bucket': BUCKET, 'sources': {}}
    index['sources'][name] = {'provenance': provenance, 'page_count': len(pages), 'objects': saved,
        'local_source_file': False, 'verification': 'full S3 readback SHA-256'}
    INDEX.write_text(json.dumps(index, indent=2)+'\n', encoding='utf-8', newline='\n')
    print(name, len(pages), 'pages', len(pdf), 'PDF bytes; S3 readback verified', flush=True)


def main():
    sys.stdout.reconfigure(encoding='utf-8')
    parser = argparse.ArgumentParser()
    parser.add_argument('action', choices=['report', 'url', 'show', 'reextract'])
    parser.add_argument('name')
    parser.add_argument('--url')
    parser.add_argument('--pages', default='1:999')
    args = parser.parse_args()
    session = boto3.Session(profile_name='quantyra')
    assert session.client('sts').get_caller_identity()['Account'] == OWNER
    s3 = session.client('s3', region_name='us-east-1')
    if args.action == 'report':
        key = 'uap-release-analysis/official-releases/release_06/documents_release_06_sept_18_2026.zip'
        # Frozen acquisition version; never silently switch to a newer bundle.
        bundle_rows = json.loads((ROOT.parent.parent / 'uap-release-analysis/docs/research/ufo-releases-05-07-s3-receipts.json').read_text())['objects']
        receipt = next(r for r in bundle_rows if r['s3_uri'].endswith('/'+key))
        version = receipt['version_id']
        remote = S3Range(s3, key, version, receipt['bytes'])
        with zipfile.ZipFile(remote) as z:
            member = next(n for n in z.namelist() if 'DOW-UAP-D153' in n and '__MACOSX' not in n)
            pdf = z.read(member)  # ZIP decompression and CRC verified in memory.
        archive(s3, args.name, pdf, {'bundle_key': key, 'bundle_version': version,
            'bundle_sha256': receipt['sha256'], 'member': member})
    elif args.action == 'url':
        js = "(async()=>{let r=await fetch(process.argv[1]);if(!r.ok)throw Error('HTTP '+r.status);process.stdout.write(Buffer.from(await r.arrayBuffer()))})().catch(e=>{console.error(String(e));process.exitCode=1})"
        result = subprocess.run(['node', '-e', js, args.url], capture_output=True, check=True)
        archive(s3, args.name, result.stdout, {'url': args.url})
    elif args.action == 'reextract':
        index = json.loads(INDEX.read_text())
        source = index['sources'][args.name]['objects']['pdf']
        data = s3.get_object(Bucket=BUCKET, Key=source['key'], VersionId=source['version_id'],
            ExpectedBucketOwner=OWNER)['Body'].read()
        assert hashlib.sha256(data).hexdigest() == source['sha256']
        pages = [p.extract_text() or '' for p in PdfReader(io.BytesIO(data)).pages]
        extracted = json.dumps(pages, ensure_ascii=False).encode('utf-8')
        key = f'{PREFIX}/{args.name}.pages.json'
        result = s3.put_object(Bucket=BUCKET, Key=key, Body=extracted, ContentType='application/json',
            ServerSideEncryption='AES256', ExpectedBucketOwner=OWNER)
        version = result['VersionId']
        check = s3.get_object(Bucket=BUCKET, Key=key, VersionId=version,
            ExpectedBucketOwner=OWNER)['Body'].read()
        digest = hashlib.sha256(extracted).hexdigest()
        assert hashlib.sha256(check).hexdigest() == digest
        entry = index['sources'][args.name]
        entry['objects']['pages.json'] = dict(key=key, version_id=version, bytes=len(extracted), sha256=digest)
        entry['extraction'] = 'pypdf 6.21.2 with fonttools 4.67.0; PDF source unchanged'
        INDEX.write_text(json.dumps(index, indent=2)+'\n', encoding='utf-8', newline='\n')
        print(args.name, 'font-aware text re-extraction verified', flush=True)
    else:
        source = json.loads(INDEX.read_text())['sources'][args.name]['objects']['pages.json']
        data = s3.get_object(Bucket=BUCKET, Key=source['key'], VersionId=source['version_id'],
            ExpectedBucketOwner=OWNER)['Body'].read()
        assert hashlib.sha256(data).hexdigest() == source['sha256']
        pages = json.loads(data)
        lo, hi = map(int, args.pages.split(':'))
        for i in range(lo-1, min(hi, len(pages))):
            print(f'\n[PDF page {i+1}]\n{pages[i]}')


if __name__ == '__main__':
    main()
