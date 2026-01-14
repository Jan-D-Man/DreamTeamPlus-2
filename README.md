README -- Estructura del projecte

# Simulació orbital i producció fotovoltaica

Aquest projecte modelitza l’òrbita terrestre mitjançant mètodes numèrics i utilitza aquesta informació per calcular la radiació solar incident i la producció energètica d’una placa fotovoltaica al llarg de l’any.

## Estructura del projecte


Entrega2plaquessolars/

    Definicions.py

    ResolucioEDO/
        Euler.py
        RK2.py
        RK4.py
        Error.py
        3solucionsnum.py

    ResolucioProblema/
        TrajectoriaSolar.py
        PotDiariaAngles.py
        PotAnualAnglesAzimut.py
        PotAnualAnglesBeta.py
        PotModeGirasol.py
        PicsPotenciaanual.py


## Descripció dels fitxers

### `Definicions.py`

Conté les constants físiques, paràmetres orbitals i dades del sistema fotovoltaic utilitzades a tot el projecte.

### Carpeta `ResolucioEDO/`

Inclou els mètodes numèrics per resoldre l’equació diferencial de l’òrbita:

* `Euler.py`: mètode d’Euler.
* `RK2.py`: mètode de Runge–Kutta d’ordre 2.
* `RK4.py`: mètode de Runge–Kutta d’ordre 4 (mètode principal).
* `Error.py`: anàlisi de l’error numèric dels mètodes.
* `3solucionsnum.py`: comparació de les solucions obtingudes amb cada mètode.

### Carpeta `ResolucioProblema/`

Scripts per al càlcul de la posició solar i la producció energètica:

* `TrajectoriaSolar.py`: trajectòria aparent del Sol al llarg de l’any.
* `PotDiariaAngles.py`: potència diària segons orientació i inclinació.
* `PotAnualAnglesAzimut.py`: energia anual en funció de l’azimut.
* `PotAnualAnglesBeta.py`: energia anual en funció de la inclinació.
* `PotModeGirasol.py`: producció amb seguiment solar.
* `PicsPotenciaanual.py`: anàlisi dels pics de potència anuals.

