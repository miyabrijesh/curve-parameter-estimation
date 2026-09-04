# Parametric Curve Parameter Estimation

This repository estimates the unknown parameters in the assigned curve:

\[
x=t\cos(\theta)-e^{M|t|}\sin(0.3t)\sin(\theta)+X
\]

\[
y=42+t\sin(\theta)+e^{M|t|}\sin(0.3t)\cos(\theta)
\]

The CSV rows are unordered, so the method does not assume a row number is a value of \(t\).

## Final values

| Parameter | Estimated value |
| --- | ---: |
| \(\theta\) | \(0.523598303\) radians (\(29.999973^\circ\)) |
| \(M\) | \(0.029999997\) |
| \(X\) | \(54.999998\) |

Rounded answer:

\[
\boxed{\theta=\pi/6,\quad M=0.03,\quad X=55}
\]

### Desmos-ready expression

```text
(t*cos(0.523598776)-e^(0.03*abs(t))*sin(0.3t)*sin(0.523598776)+55,
42+t*sin(0.523598776)+e^(0.03*abs(t))*sin(0.3t)*cos(0.523598776))
```

Use \(6<t<60\).

## Method

For a candidate \(\theta\) and \(X\), undo its rotation and translation for every observed point:

\[
t_i=\cos(\theta)(x_i-X)+\sin(\theta)(y_i-42)
\]

\[
b_i=-\sin(\theta)(x_i-X)+\cos(\theta)(y_i-42)
\]

The curve requires \(b_i=e^{Mt_i}\sin(0.3t_i)\). The fit therefore minimizes:

\[
r_i=b_i-e^{Mt_i}\sin(0.3t_i)
\]

The solver uses differential evolution for a global bounded search, then bounded nonlinear least-squares to refine it. This avoids a lucky initial guess and works directly with unordered samples.

## Verification

- 1,500 supplied data points
- Recovered \(t\): 6.0494 to 59.9952, within the required range
- Mean absolute residual: \(2.56\times10^{-6}\)
- RMSE: \(3.49\times10^{-6}\)

Those errors are at CSV rounding level, strongly supporting the rounded answer above.

## Run

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
python estimate_parameters.py --data xy_data.csv
```

This generates `output/fit_summary.json` and `output/curve_fit.png`.
