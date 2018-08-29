#
# Copyright (C) 2018  Leo Singer
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 2 of the License, or (at your
# option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
# Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
#
"""Make common LAL datatypes picklable."""
import copyreg

import lal


numpy_to_lal_types = {'char': 'CHAR',
                      'int16': 'INT2',
                      'int32': 'INT4',
                      'int64': 'INT8',
                      'uint16': 'UINT2',
                      'uint32': 'UINT4',
                      'uint64': 'UINT8',
                      'float32': 'REAL4',
                      'float64': 'REAL8',
                      'complex64': 'COMPLEX8',
                      'complex128': 'COMPLEX16'}


def pickle_gps(obj):
    return lal.LIGOTimeGPS, (obj.gpsSeconds, obj.gpsNanoSeconds)


def pickle_unit(obj):
    return lal.Unit, (str(obj),)


def vector(data):
    lal_type = numpy_to_lal_types[data.dtype.name]
    creator = getattr(lal, 'Create{}Vector'.format(lal_type))
    result = creator(len(data))
    result.data = data
    return result


def pickle_vector(obj):
    return vector, (obj.data,)


def series(attrs):
    lal_type = numpy_to_lal_types[attrs['data'].data.dtype.name]
    kind = 'Frequency' if 'deltaF' in attrs else 'Time'
    creator = getattr(lal, '{}{}Series'.format(lal_type, kind))
    result = creator()
    for key, value in attrs.items():
        setattr(result, key, value)
    return result


def pickle_series(obj):
    attrs = {'name': obj.name, 'epoch': obj.epoch, 'f0': obj.f0,
             'sampleUnits': obj.sampleUnits, 'data': obj.data}
    if hasattr(obj, 'deltaF'):
        attrs['deltaF'] = obj.deltaF
    else:
        attrs['deltaT'] = obj.deltaT
    return series, (attrs,)


copyreg.pickle(lal.LIGOTimeGPS, pickle_gps)
copyreg.pickle(lal.Unit, pickle_unit)
for datatype in numpy_to_lal_types.values():
    clazz = getattr(lal, '{}Vector'.format(datatype), None)
    if clazz:
        copyreg.pickle(clazz, pickle_vector)
    clazz = getattr(lal, '{}FrequencySeries'.format(datatype), None)
    if clazz:
        copyreg.pickle(clazz, pickle_series)
    clazz = getattr(lal, '{}TimeSeries'.format(datatype), None)
    if clazz:
        copyreg.pickle(clazz, pickle_series)
