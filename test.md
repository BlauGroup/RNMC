---
title: TESTING
layout: template
filename: test
--- 

there are 3 tables in the initial state database:
```
CREATE TABLE initial_state (
    site_id            INTEGER NOT NULL PRIMARY KEY,
    degree_of_freedom  INTEGER NOT NULL
);
```

```
CREATE TABLE trajectories (
    seed               INTEGER NOT NULL,
    step               INTEGER NOT NULL,
    time               REAL NOT NULL,
    site_id_1          INTEGER NOT NULL,
    site_id_2          INTEGER NOT NULL,
    interaction_id     INTEGER NOT NULL
);
```


```
CREATE TABLE factors (
    one_site_interaction_factor      REAL NOT NULL,
    two_site_interaction_factor      REAL NOT NULL,
    interaction_radius_bound         REAL NOT NULL,
    distance_factor_type             TEXT NOT NULL
);
```

`distance_factor_type` specifies how to compute interaction propensities for two site interactions as a function of distance. Currently the accepted values are `linear` and `inverse_cubic`.