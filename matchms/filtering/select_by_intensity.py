import numpy


def select_by_intensity(spectrum_in, intensity_from=10.0, intensity_to=200.0):

    if spectrum_in is None:
        return None

    spectrum = spectrum_in.clone()

    assert intensity_from <= intensity_to, "'intensity_from' should be smaller than or equal to 'intensity_to'."

    condition = numpy.logical_and(intensity_from <= spectrum.intensities, spectrum.intensities <= intensity_to)

    spectrum.mz = spectrum.mz[condition]
    spectrum.intensities = spectrum.intensities[condition]

    return spectrum
