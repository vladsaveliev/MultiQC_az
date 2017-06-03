#!/usr/bin/env python
""" MultiQC command line options - we tie into the MultiQC
core here and add some new command line parameters. """

import click

az_metadata = click.option('--az-metadata', 'az_metadata',
    type = str,
    help = 'File with AZ specific reporting metadata'
)
