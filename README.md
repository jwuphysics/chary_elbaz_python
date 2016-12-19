# A python version of the Chary & Elbaz (2001) code

All files in the directory `chary_elbaz` are housed [here](http://david.elbaz3.free.fr/astro_codes/chary_elbaz.html) [(Chary & Elbaz 2001, ApJ, 556, 562)](http://adsabs.harvard.edu/abs/2001ApJ...556..562C). Please cite their paper if you use this program in your research.

This code is a simple python port of their IDL script for converting far-infrared flux density to total infrared luminosity (8-1000 micron). I use the `simps` integration method from SciPy, which slightly differs from the `INT_TABULATED` integration routine in IDL. In my own usage, I've found differences of about 0.03% between the two.

An example of how to use this code:

```python
import chary_elbaz as ce

# find infrared luminosity of a 3.5 mJy source at 160 microns, at z=1
L_IR = ce.quick_fit_SED(z=1, S_nu=3500, wavelength=160)

# or equivalently
import astropy.units as u
L_IR = ce.fit_SED(z=1, S_nu=3500*u.uJy, wavelength=160*u.um)

print(L_IR * 1e-10) 
# 16.1544513536
```
We see that the fit infrared luminosity is 16.2e10 solar luminosities. Note that the returned quantity has no explicit units.
