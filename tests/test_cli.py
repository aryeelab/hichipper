import pytest
from click.testing import CliRunner
from hichipper import cli
import md5


def file_checksums_equal(file1, file2):
    with open(file1) as f:
        checksum1 = md5.new(f.read()).digest()
    with open(file2) as f:
        checksum2 = md5.new(f.read()).digest()
    return checksum1==checksum2 


def test_loops_output():
	runner = CliRunner()
	result = runner.invoke(cli.main, ['--out', 'output1', 'example.yaml'])
	assert file_checksums_equal('correct_output/co.intra.loop_counts.bedpe', 'output1/test_sample1.intra.loop_counts.bedpe')
	assert file_checksums_equal('correct_output/co.intra.loop_counts.bedpe', 'output1/test_sample2.intra.loop_counts.bedpe')