# procomom multiserver regression fixtures

Four small replicated-station models built to cover the `mi>1` path in
`procomom`/`promom`/`camarg`. They live here rather than in `models/` so the
main corpus, which every solver sweep iterates over, stays the size it was.

| file | shape | why |
|---|---|---|
| `m1.qn` | `mi=3` at the first station, R=2 | replicated station away from the reference |
| `m2.qn` | `mi=4` at the LAST station, R=2 | with the reference rotated to position M, this exercises the replicated-REFERENCE path for every station in turn -- the defect that needed the `M+1` stride |
| `m3.qn` | every station replicated, R=2 | multiplicities interacting across the whole PC row |
| `m4.qn` | R=3, `Z != 0`, mixed `mi` | think times plus multiplicities together |

All four agree exactly across the three independent routes:

```
./bin/ca <f> -q        (summed over classes)
./bin/camarg <f> -q    == ca, and -P bit-for-bit == procomom -P
./bin/procomom <f> -q  == ca
```

They exist because the real replicated models in `models/` are impractical as
a check here: `10_diverse` needs about an hour (basis 1320, 560 steps) and
`11_swapped` is degenerate, so `procomom` perturbs and its coefficients
explode. See `procomom/README.md` for both cost regimes.
