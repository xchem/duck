import os
import tempfile
from urllib.parse import urlparse
from urllib.request import pathname2url
import boto3
from botocore.exceptions import ClientError


def get_bucket_and_key_from_url(url):
    if not url.startswith("s3://"):
        raise ValueError(
            "cannot parse S3 URL: only the s3:// scheme is supported.")
    result = urlparse(url)
    return result.netloc, result.path.lstrip('/')


def get_file_name_from_url(url):
    if not url.startswith("s3://"):
        raise ValueError(
            "cannot parse S3 URL: only the s3:// scheme is supported.")
    _, fname = url.rsplit('/', 1)
    return fname


def upload_file_to_s3(file, url):
    bucket_name, key = get_bucket_and_key_from_url(url)
    if not check_if_bucket_exists(bucket_name):
        raise ValueError(
            "cannot upload file: S3 bucket "+bucket_name+" does not exist.")
    s3 = boto3.resource("s3")
    bucket = s3.Bucket(bucket_name)
    bucket.upload_file(file, key)


def download_file_from_s3(url, path):
    bucket_name, key = get_bucket_and_key_from_url(url)
    s3 = boto3.resource("s3")
    bucket = s3.Bucket(bucket_name)
    try:
        bucket.download_file(key, path)
    except ClientError as err:
        err_code = err.response["Error"]["Code"]
        if err_code == "404":
            raise ValueError(
                "failed to download file from S3: {url} does not exist.")
        elif err_code == "403":
            raise ValueError(
                "failed to download file from S3: access denied.")
        else:
            raise ValueError("failed to download file from S3: "+err_code)


def download_s3_to_temp_file(
        url, prefix=None, dir= None):
    ext = '.' + os.path.basename(url).split('.')[0]
    _, temp_path = tempfile.mkstemp(suffix=ext, prefix=prefix, dir=dir)
    try:
        download_file_from_s3(url, temp_path)
    except ValueError as err:
        os.remove(temp_path)
        raise ValueError(str(err))
    return temp_path

def check_if_bucket_exists(bucket):
    s3 = boto3.resource("s3")
    try:
        s3.meta.client.head_bucket(Bucket=bucket)
    except ClientError as err:
        error_code = int(err.response["Error"]["Code"])

        if error_code == 404:
            return False
    return True

def write_data_to_s3(data, url):
    bucket, key = get_bucket_and_key_from_url(url)
    if not check_if_bucket_exists(bucket):
        raise ValueError(
            "cannot upload file: S3 bucket "+bucket+" does not exist.")
    key = pathname2url(key)
    s3 = boto3.resource("s3")
    bucket = s3.Bucket(bucket)
    bucket.put_object(Key=key, Body=data)


def copy_directory_to_s3(directory_path, url):
    if not url.endswith('/'):
        raise ValueError(
            "cannot upload files: S3 bucket {url} key has incorrect format.")
    # Iterate over the files in the directory
    files = os.listdir(directory_path)
    num_files = len(files)
    for filename in files:
        local_path = os.path.join(directory_path,filename)
        if os.path.isdir(local_path):
            continue
        remote_path = os.path.join(url,filename)
        data = open(local_path, 'rb')
        write_data_to_s3(data, remote_path)





