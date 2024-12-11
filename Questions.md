# Questions du TP Poisson

Réponses aux questions du TP poisson.

## Exercice 1 

## Exercice 2

## Exercice 3

### Question 1
Pour utiliser `BLAS` et `LAPACK`, nous devons déclarer et allouer les matrices sous forme d'**un tableau 1D contigu** (étant considéré comme un 2D).
Nous utiliserons pour cela `LAPACK_ROW_MAJOR` en C pour déclarer les tableaus en `Ligne-Major`.

Exemple:
```c++
double* A = (double*)malloc(m * n * sizeof(double));
```

### Question 2

La constante `LAPACK_COL_MAJOR` fait référence à **l'ordre de stockage en mémoire**.
Ici les matrices seront stockées en `Colonne-Major` (plutôt pour Fortran que C).

### Question 3

La dimension principale (`ld`, ou `lda` pour leading dimension of A) permet d'indiquer comment une matrice est stockée en mémoire.
Elle permet d'indiquer le **nombre d'éléments entre chaque ligne** (en colonne-major) ou chaque colonne (en ligne-major).


