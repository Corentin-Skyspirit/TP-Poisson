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

### Question 4

La fonction `dgbmv` **effectue une multiplication matrice-vecteur** pour une matrice bande.

### Question 5

La fonction `dgbrtf` effectue la **factorisation LU** d'une matrice bande avec **pivotement partiel**.

### Question 6

La fonction `dgbrts` **résout un système d'équation linéaires** en utilisant la factorisation LU d'une matrice bande préalablement calculée par `dgbrtf`.

### Question 7

La fonction `dgbsv` effectue à **la fois la factorisation LU et la résolution d'un système linéaire**.
Elle permet de faire en une fonction `dgbrtf` et `dgbrts`.

### Question 8

Soit R le résidu : $R = b - A\^x$

L'erreur arrière est : 
$$\frac{||b - A\^x||}{||A||||\^x||}$$

L'erreur avant est : 
$$\frac{||k||}{||b||} = \frac{||b - A\^x||}{||b||} = \frac{||x - \^x||}{||x||}$$

Pour calculer la norme du vecteur : $n = sqrt(ddot(\&n, x, 1, x, 1))$

Pour calculer $||x - \^x||$ : $daxpy(\&n, x, 1, x, 1) = \^x-x = x$

et $nx = norme(x)$ avec $res = n / nx$

## Exercice 4

### Question 2 & 3

Afin de valider les valeurs obtenues dans la matrice, il faut multiplier un **vecteur unitaire avec AB** avec la méthode `dgbmv`.
Il suffira ensuite de vérifier que le résultat correspond à un **vecteur nul avec un 1 à chaque extrémités**.

## Exercice 5

### Question 2

Les performances de `dgbsv` sont nettement supérieures aux performances de `dgbtrs`. Cela semble logique au vu de la complexité de `dgbsv` qui est inféireure à `dgbtrs` car elle s'effectue sur une martrice bande.
