#!/usr/bin/env python

"""
-------------------------------------------------------------------------------------------------------------------------------------------
MIT License

Copyright (c) 2024 George P. Tiley

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the “Software”), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

Contact: george.tiley@gmail.com
-------------------------------------------------------------------------------------------------------------------------------------------

estploidy uses a command-line interface to interact with the estploidy package. A successful installation should allow you to access all functions and their help with:

estploidy --help
"""

import click
import time
from datetime import datetime
import logging

#--------------------------------#
# CLI ENTRY POINT
#--------------------------------#

@click.group()
def cli():
    current_time = datetime.now()
    time_string = current_time.strftime('%Y-%m-%d-%H-%M-%S')
    logfile = time_string + 'estploidy.log'
    print(f'Start estploidy at {current_time}')
    logging.basicConfig(
       filename = logfile,
       format='%(asctime)s - %(levelname)s - %(message)s',
       filemode = 'w' 
    )
    print(f'Created logfile {logfile}')



@cli.command(context_settings={'help_option_names': ['-h','--help']})

@click.argument('sample_sheet',type=str)
@click.option('-v', '--vcf_file', type=str, default='dummy.vcf', required=True,
              help = 'name of the input vcf file. should not be compressed.'
)
@click.option('-i', '--imputation_method', type=str, default='drop', required=False,
              help = 'decide how to impute missing data if at all. options are: drop, mean, and popmean'
)
@click.option('-d', '--minimum_depth', type=int, default=10, required=False,
              help = 'The minimum depth of a site to be treated as data'
)
@click.option('-c', '--minimum_count', type=int, default=3, required=False,
              help = 'The minimum count of the minor allele for a site to be treated as data'
)
@click.option('-q', '--minimum_quality', type=int, default=20, required=False,
              help = 'The minimum phred-scaled genotype quality score'
)
@click.option('-o', '--output_dir', type=str, default='dummy', required=False,
              help = 'name of the directory where . will be a matrix of allele frequencies'
)

def individual_frequencies(sample_sheet, vcf_file, minimum_depth, minimum_count, minimum_quality, imputation_method, output_dir):
    from estploidy.utils import map_individuals
    from estploidy.utils import check_dir
    from estploidy.calculate_frequencies.calculate_frequencies import get_ind_freqs

    start_time = time.process_time()
    logging.info(f'Begin at {start_time}')
    logging.info(f'Checking all individuals in {sample_sheet} are present in {vcf_file}')
    ind_map = map_individuals(sample_sheet, vcf_file)
    logging.info(f'Calculating individual allele frequencies from {vcf_file}')
    if (imputation_method == 'drop'):
        if (output_dir != 'dummy'):
            check_dir(output_dir)
            logging.info(f'Matrix of allele frequencies for each individual will be written to: {output_dir}')
        get_ind_freqs(ind_map, vcf_file, minimum_depth, minimum_count, minimum_quality, output_dir)
        
    else:
        click.echo(f'Warning: Imputation method {imputation_method} is not supported. Skipping allele frequencies.')
    end_time = time.process_time()
    logging.info(f'End at {end_time}')
    compute_time = (end_time - start_time) / 60
    logging.info(f'Total compute time was {compute_time} minutes')

##----------------
## Calculate population-level allele frequencies for fst and genotype-environment association analyses
##----------------
#@click.argument('sample_sheet',type=str)
#@click.option('-v', '--vcf_file', type=str, default='dummy.vcf',
#              help = 'name of the input vcf file. should not be compressed.'
#)
#@click.option('-o', '--output_file', type=str, default='samples.vcf',
#              help = 'name of the output text file. will be a matrix of allele frequencies'
#)
#def population_frequencies(sample_sheet, vcf_file, output_file):
#    pass
