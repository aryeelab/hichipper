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
	result = runner.invoke(cli.main, ['--out', 'output1', '--peak-pad', '1000', '--skip-resfrag-pad', '--skip-qc', '--skip-diffloop', 'yaml/one.yaml'])
	assert file_checksums_equal('correct_output/goal.loop_counts.bedpe', 'output1/d1.filt.intra.loop_counts.bedpe')
