# Sky Model

A sky model consists of one of more named sources. Each source consists of one or more emission components with properties that describe their respective position on the sky, type, and emission spectrum. This information can be used to describe both bright calibrator sources (like Cygnus A, for example) or an entire field of view of fainter, smaller, point sources. And either can be used as a model to obtain a calibration solution.

This page describes the format used in sky model files, passed to calibrate using the `--model` parameter. This page also describes how some of these sky model parameters are interpreted and processed.

```@contents
Pages = ["skymodel.md"]
Depth = 5
```

## skymodel fileformat 1.1

### Syntax

Keywords and tokens are required to be separated by white space (one more or of spaces, tabs, new lines), with the exception of braces (`{` or `}`) which are understood as complete keywords.

Single quotes (`'`) or double quotes (`"`) can be used around a strings containing white space, for example in source names like `"Crab Nebula"` or `'Fornax A'`, so long as the opening and closing quotes match. There are no escape characters for within quoted strings.

The format is otherwise not white-space sensitive: new lines, spaces and tabs can be used to format the file as desired.

Comments are indicated by a `#` symbol and will exclude all proceeding text until the next new line from further processing.

!!! note
    The function `gettoken()` can be referenced for further clarity about token parsing.

```@docs
gettoken
```

### Format

The sky model file must include the following header as the first non-whitespace text:

    skymodel fileformat 1.1

Following this, any number of sources may be provided:

    source {
        name "Source name"
        component { ... } [component { ... }, ...]
    }

A **source** provides a representation of a radio object. For a typical point-like source on the sky, a source will have just one component. For larger or more complex sources, a source may have hundreds or even thousands of components.

Each source must have a `name`, and optionally 0 or more components.

    component {
        position 00h00m00.0s 00d00m00.0s
        type {point | gaussian}
        [shape 150 120 30]
        [measurement {...}, ...] | sed{...}]
    }

A **component** describes emission originating from a single direction in the sky. Each component must have a position, a type and spectrum.

The `position` describes the J2000 RA and Dec coordinates of the emission and must be in the sexagesimal format illustrated above, using the 'h', 'm' and 's' characters for RA and 'd', 'm', 's' for Dec. No other format is accepted.

The `type` of the component can be a point source or a Gaussian. A point source is represented analytically in the measurement equation as a Dirac delta function.

A gaussian source requires the additional parameter `shape` which describes the major and minor axis (arcseconds) and a position angle (degrees). _**ToDo:** Which direction is the angle measured from?_

Each component must also contain spectrum information. The spectrum describes the flux of the source as a function of frequency. This is provided in one of two ways: either as one or more measurements at specified frequencies, or as a spectral energy distrbution (SED).

A **measurement** is of the format:

    measurement {
      frequency 100 MHz
      fluxdensity Jy 100 2 2 0
    }

The `frequency` field describes the frequency at which the meausurement was taken. The unit `MHz` must be appended, and this is the only unit currently accepted. Other units like `GHz` will throw an error.

`fluxdensity` describes the measured flux of the component at the specified frequency in units of Jansky in Stokes I, Q, U, V, respectively. The unit `Jy` must be provided, and is the only unit currently accepted.

!!! note
    The flux density is the total/integrated flux of the source [Jy], not the peak flux [Jy / beam]. For point sources these values are identical, but for Gaussian sources they are not.

!!! note
    [See below](#measurement-interpolation-rules) for the rules for how flux density is interpolated in frequency when using measurements.

An **SED** is an alternative to providing one more more measurements. It allows for complete control over how the spectrum is interpolated in frequency by specifying arbitrary Taylor terms describing the curve:

```math
\log(S) = \log(S_0) + \alpha \log(\nu - \nu_0) + \beta \log(\nu - \nu_0)^2 + ...
```

The format for specifying an SED is as follows:

    sed {
      frequency 154 MHz
      fluxdensity Jy 100 0 0 0
      spectral-index { -0.7 [...] }
    }

The `frequency` field specifies the reference frequency ``\nu_0``. Only `MHz` is accepted.

The `fluxdensity` specifies the flux density ``S_0`` measured at the reference frequency in units of Jansky in Stokes I, Q, U, V, respectively. The unit `Jy` must be provided, and is the only unit currently accepted. The flux density is the total/integrated flux of the source.

The spectral index describes the Taylor terms of the curve in log/log space, i.e. ``\alpha``, ``\beta``, etc. At least one Taylor term must be provided, however there is no limit to the number of coefficients that are accepted.

!!! note
    The SED curve applies to all of the Stokes parameters, and there is not currently any way to specify unique curves for each.

## Example skymodels

### A single point source with measurements

    skymodel fileformat 1.1
    source {
      name "J035857+102702"
      component {
        type point
        position 03h58m57.7099s +10d27m17.892s
        measurement {
          frequency 80 MHz
          fluxdensity Jy 46.75907 0 0 0
        }
        measurement {
          frequency 100 MHz
          fluxdensity Jy 39.73867 0 0 0
        }
        measurement {
          frequency 120 MHz
          fluxdensity Jy 34.80545 0 0 0
        }
      }
    }

### A single point source with SED

    source {
      name "Crab"
      component {
        type point
        position 05h34m28.1s 22d02m09s
        sed {
          frequency 154 MHz
          fluxdensity Jy 87.0601723 0 0 0
          spectral-index { -0.22 0.00 }
        }
      }
    }

### A multicomponent source with mixed component types

    source {
      name "J232013-132102A"
      component {
        type gaussian
        position 23h20m10.1296s -13d19m47.316s
        shape 185.821633 159.627594 -89.000000
        sed {
          frequency 154 MHz
          fluxdensity Jy 87.0601723 0 0 0
          spectral-index { -0.22 0.00 }
        }
      }
      component {
        type point
        position 23h20m17.9149s -13d22m31.836s
        measurement {
          frequency 80 MHz
          fluxdensity Jy 0.27617 0 0 0
        }
        measurement {
          frequency 240 MHz
          fluxdensity Jy 0.15225 0 0 0
        }
      }
    }

## Measurement interpolation rules

When specifying a component spectrum using measurements, it is usually necessary for calibrate to interpolate in frequency between the specified measurements. These rules follow exactly the rules used by the original claibrate, as well as handling all edge cases.

The rules for how this works are best described by the documentation of the implementing function:

```@docs
stokes
```
