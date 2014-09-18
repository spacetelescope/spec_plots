Edge trimming checks include:
1.  DQ flag > 0.
2.  Flux is > 5 * median flux (median flux calculation ignores any fluxes = 0.0).
3.  Flux error is consistent with zero flux to within 1-sigma.
4.  Flux value is exactly 0.0 (STIS) or < 0.0 (COS).

Two different sets of start/end indexes are returned.  The first includes the DQ test, the second does not.  Some spectra have all fluxes with DQ > 0, in which case it falls back to the test that does not check DQ flags to determine the start/end indexes of the spectrum.

DQ flag > 0 is way too liberal of a cut, need to be more selective on what DQ flags I remove, at least for STIS.
