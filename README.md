# Hill Cipher — MPI C++ Implementation

Parallel/distributed implementation of the Hill cipher written in C++ using MPI.

This repository contains `hill_cipher_mpi.cpp`, a project that demonstrates how to
encrypt and decrypt text with the Hill cipher while distributing work across
multiple processes with MPI for performance and learning purposes.

## Features

- Matrix-based Hill cipher (encryption/decryption)
- Parallel workload distribution using MPI (message passing)
- Example build and run instructions for Linux/macOS and Windows

## Prerequisites

- A C++ compiler with MPI support (e.g., `mpic++` from Open MPI or MPICH)
- MPI runner: `mpirun` / `mpiexec`
- Basic knowledge of compiling C++ and running MPI programs

On Windows you can use MS-MPI or an MPI distribution that supplies `mpic++` and
`mpiexec`.

## Build

From the repository root (where `hill_cipher_mpi.cpp` is located) run:

```bash
mpic++ -O2 -std=c++17 -o hill_cipher_mpi hill_cipher_mpi.cpp
```

On Windows (PowerShell or CMD) with an MPI-enabled `mpic++` installed the same
command should work and produce `hill_cipher_mpi.exe`.

## Run

Replace `<num_procs>` and the program arguments below with values your program expects.

Linux / macOS:
```bash
mpirun -np <num_procs> ./hill_cipher_mpi [args]
```

Windows (MS-MPI or mpiexec):
```powershell
mpiexec -n <num_procs> hill_cipher_mpi.exe [args]
```

Common example patterns (adjust to your program's CLI):

- Encrypt a file: `mpirun -np 4 ./hill_cipher_mpi encrypt key_matrix.txt plaintext.txt ciphertext.txt`
- Decrypt a file: `mpirun -np 4 ./hill_cipher_mpi decrypt key_matrix.txt ciphertext.txt decrypted.txt`

If `hill_cipher_mpi.cpp` uses interactive input or different command-line flags,
update the commands above to match its interface.

## Usage / Notes

- The Hill cipher operates on blocks sized to the key matrix (e.g., 2x2, 3x3).
- Ensure the key matrix is invertible modulo the alphabet size (commonly 26).
- The program parallelizes either by splitting the input text into chunks or by
  distributing matrix operations; check the source for the exact approach.

## Key Files

- `hill_cipher_mpi.cpp` — main implementation (encryption/decryption + MPI logic)

## Example

1. Build:

```bash
mpic++ -O2 -std=c++17 -o hill_cipher_mpi hill_cipher_mpi.cpp
```

2. Run with 4 processes (example):

```bash
mpirun -np 4 ./hill_cipher_mpi encrypt key.txt plain.txt cipher.txt
```

## Troubleshooting

- If you see MPI errors, verify your MPI installation and that `mpic++` matches
  the `mpirun`/`mpiexec` distribution.
- For Windows, prefer launching from the Developer Command Prompt or ensure
  environment variables point to the MPI bin folder.

## License & Credits

Add a license file if you wish to open-source this code (e.g., MIT, Apache-2.0).

---

If you'd like, I can:

- Adjust the README to match the exact CLI and behavior of `hill_cipher_mpi.cpp` (I can
  inspect the source and fill in precise examples), or
- Add a sample key and input files and a short test harness.

Tell me which option you prefer and I'll update the README accordingly.
"# Distributed-Hill-Cipher-C-MPI-" 
