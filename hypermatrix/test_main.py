#!/usr/bin/env python3
# File: test_main.py

import unittest
import os
import sys
import json
import shutil
from unittest.mock import patch
from io import StringIO
import logging

# Import the main module
import main

class TestMain(unittest.TestCase):
    def setUp(self):
        # Define paths to your actual test data
        self.test_data_dir = os.path.expanduser('~/test_data')
        self.output_dir = os.path.join(self.test_data_dir, 'output')
        os.makedirs(self.output_dir, exist_ok=True)

        # Paths to your data files
        self.methy_dir = self.test_data_dir  # Assuming methylation data is here
        self.hic_file = os.path.join(self.test_data_dir, 'GSE158007_DKO1_Hi-C_allValidPairs.hic')
        self.genome_fa = os.path.join(self.test_data_dir, 'hg38.fa.gz')
        self.chrom_sizes = os.path.join(self.test_data_dir, 'hg38.chrom.sizes')

        # Set up test config
        self.config_path = os.path.join(self.output_dir, 'config.json')
        self.original_config_path = main.config_path  # Save original config path
        main.config_path = self.config_path  # Redirect to test config

        # Create a sample config using your actual paths
        sample_config = {
            "bam_directory": self.hic_file,
            "methy_directory": self.methy_dir,
            "output_directory": self.output_dir,
            "reference_genome": "hg38",
            "genome_fa": self.genome_fa,
            "chrom_file": self.chrom_sizes,
            "resolutions": ["40000:40kb"],
            "impute": False,
            "cluster_compartments": False,
            "cumulant": False,
            "iterations": 400,
            "chromosomes": [str(i) for i in range(1, 23)],
            "mappability_threshold": 0.5,
            "normalization": "NONE",
            "data_type": "oe",
            "correlation": True,
            "hic_cell_type1_url": "",
            "hic_cell_type2_url": ""
        }

        with open(self.config_path, 'w') as f:
            json.dump(sample_config, f, indent=4)

        # Set up logging to capture outputs for testing
        self.log_stream = StringIO()
        logger = logging.getLogger()
        logger.handlers = []  # Clear existing handlers
        handler = logging.StreamHandler(self.log_stream)
        logger.addHandler(handler)
        logger.setLevel(logging.INFO)

    def tearDown(self):
        # Clean up and restore original config path
        main.config_path = self.original_config_path
        logging.shutdown()

    @patch('main.run_shell_script')
    def test_abcluster_with_args(self, mock_run_script):
        # Mock sys.argv to simulate command-line arguments
        with patch('sys.argv', ['main.py', 'ABcluster',
                                '-m', self.methy_dir,
                                '-i', self.hic_file,
                                '-o', self.output_dir,
                                '-g', 'hg38',
                                '-r', '40000:40kb',
                                '--genome_fa', self.genome_fa,
                                '--chrom_sizes', self.chrom_sizes]):
            main.main()
            # Capture the log output
            log_output = self.log_stream.getvalue()
            # Check that the ABcluster module ran
            self.assertIn("Running ABcluster module", log_output)
            # Check that the configuration was updated
            self.assertIn("Configuration updated", log_output)
            # Ensure the correct script was called
            mock_run_script.assert_called_once()
            script_name = mock_run_script.call_args[0][0]
            self.assertIn("standard_pipeline", script_name)

    @patch('main.run_shell_script')
    def test_preprocess_nomehic(self, mock_run_script):
        # Mock sys.argv to simulate command-line arguments
        with patch('sys.argv', ['main.py', 'preprocess',
                                '--nomehic',
                                '--input_dir', self.methy_dir,
                                '--output_dir', self.output_dir,
                                '--ref_genome', 'hg38',
                                '--genome_fa', self.genome_fa,
                                '--chrom_sizes', self.chrom_sizes]):
            main.main()
            # Capture the log output
            log_output = self.log_stream.getvalue()
            # Check that the preprocess module ran
            self.assertIn("Preprocess BAM files for analysis", log_output)
            # Ensure the correct script was called
            mock_run_script.assert_called_once()
            script_name = mock_run_script.call_args[0][0]
            self.assertIn("nomehic_preprocess.sh", script_name)

    def test_load_config(self):
        # Test loading configuration
        config = main.load_config()
        self.assertIn('reference_genome', config)
        self.assertEqual(config['reference_genome'], 'hg38')

    def test_save_config(self):
        # Test saving configuration
        config = main.load_config()
        config['test_param'] = 'test_value'
        main.save_config(config)
        # Reload and verify
        config_loaded = main.load_config()
        self.assertIn('test_param', config_loaded)
        self.assertEqual(config_loaded['test_param'], 'test_value')

    @patch('sys.argv', ['main.py', 'unknown_command'])
    def test_unknown_command(self):
        # Capture stderr to check for error message
        with patch('sys.stderr', new=StringIO()) as fake_err:
            with self.assertRaises(SystemExit):
                main.main()
            output = fake_err.getvalue()
            self.assertIn("invalid choice", output)

    @patch('sys.argv', ['main.py'])
    def test_no_command(self):
        # Capture stderr to check for error message
        with patch('sys.stderr', new=StringIO()) as fake_err:
            with self.assertRaises(SystemExit):
                main.main()
            output = fake_err.getvalue()
            self.assertIn("the following arguments are required: command", output)

if __name__ == '__main__':
    unittest.main()
