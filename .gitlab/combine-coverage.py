#!/usr/bin/env python
#
# Copyright (C) 2020  Leo Singer
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
"""Combine two cobertura coverage files."""

import argparse
import lxml.etree

parser = argparse.ArgumentParser()
parser.add_argument('input1')
parser.add_argument('input2')
parser.add_argument('output')
args = parser.parse_args()

doc1 = lxml.etree.parse(args.input1)
doc2 = lxml.etree.parse(args.input2)
root1 = doc1.getroot()
root2 = doc2.getroot()
root1.attrib['lines-covered'] = str(
    int(root1.attrib['lines-covered']) +
    int(root2.attrib['lines-covered']))
root1.attrib['lines-valid'] = str(
    int(root1.attrib['lines-valid']) +
    int(root2.attrib['lines-valid']))
try:
    root1.attrib['line-rate'] = str(
        int(root1.attrib['lines-covered']) /
        int(root1.attrib['lines-valid']))
except ZeroDivisionError:
    root1.attrib['line-rate'] = '0'
root1.attrib['branches-covered'] = str(
    int(root1.attrib['branches-covered']) +
    int(root2.attrib['branches-covered']))
root1.attrib['branches-valid'] = str(
    int(root1.attrib['branches-valid']) +
    int(root2.attrib['branches-valid']))
try:
    root1.attrib['branch-rate'] = str(
        int(root1.attrib['branches-covered']) /
        int(root1.attrib['branches-valid']))
except ZeroDivisionError:
    root1.attrib['branch-rate'] = '0'
packages = root1.find('./packages')
packages.extend(root2.iterfind('./packages/package'))
doc1.write(args.output)
