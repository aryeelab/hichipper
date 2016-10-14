import pytest
from click.testing import CliRunner
from hichipper import cli

def test_preproc_run():
    runner = CliRunner()
    result = runner.invoke(cli.main, ['--out', 'output1',  'example.yaml'])
    assert not result.exception
    assert result.exit_code == 0

def test_reads_for_peak_calling():
	assert file_checksums_equal('correct_output/co.intra.loop_counts.bedpe', 'output1/test_sample1.intra.loop_counts.bedpe')
	assert file_checksums_equal('correct_output/co.intra.loop_counts.bedpe', 'output1/test_sample2.intra.loop_counts.bedpe')