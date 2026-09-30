# Tiny HAL integration fixture

This fixture is generated from tiny checked-in MAF sources rather than a committed
opaque HAL file. It is intended to quantify `halLiftover` behavior on interpretable
cases before any real bird/fish pilot.

Topology:

```text
          HAL_root
          /      \
     HAL_ancAB  HAL_ancCD
       /   \      /   \
    HAL_A HAL_B HAL_C HAL_D
```

Logical guide-tree IDs intentionally differ from HAL genome names:

```text
((A,B)ancAB,(C,D)ancCD)root;
A      -> HAL_A
B      -> HAL_B
C      -> HAL_C
D      -> HAL_D
ancAB  -> HAL_ancAB
ancCD  -> HAL_ancCD
root   -> HAL_root
```

Build the fixture locally with:

```bash
scripts/build_tiny_hal_fixture.sh
```

The script requires `maf2hal`, `halStats`, and `halLiftover` from the official
Comparative Genomics Toolkit HAL tools. CI installs these through the pinned
Bioconda package declared in `.github/workflows/ci.yml`.

The MAF sources are deliberately small and synthetic. They include blocks for a
unique collinear interval, an inverted interval, an interval absent from the
alignment, and a split/gapped interval represented by multiple shorter rows.
Genuine paralogous duplication in HAL imported from tiny hand-written MAF is not
currently treated as established; the manifest records this as an expected
limitation rather than simulating it as a real-HAL duplication.
