# Security

DiagHam is a numerical research code: it reads the files you give it and writes spectra and vectors. It does not open network connections and runs with the privileges of the user who starts it. The realistic security topics are:

- **Untrusted input files** (pseudopotential files, vector files, reference states) can crash a program (upstream reads them with fixed-size buffers in places). Do not run DiagHam on files from sources you do not trust with the rights of the account running it.
- **Build dependencies** come from your distribution or module system; this repository pins nothing that is downloaded at build time (nothing is downloaded).

To report a vulnerability — for example a way to make a program write outside its working directory — use GitHub's private vulnerability reporting on this repository ("Security" tab → "Report a vulnerability") rather than a public issue. You will get an acknowledgement; fixes are published as ordinary commits with the report credited if you wish.
