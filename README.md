# Modelos utilizados en la aplicación web *mlbiopredict*

En este repositorio se encuentra el código utilizado para desarrollar los modelos del seminario de investigación **"Desarrollo de una herramienta web para la prediccion de la actividad biologica de estructuras moleculares empleando *machine learning*"**. La página web en la que se desplegaron los modelos se encuentra disponible en este [enlace](https://mlbiopredict.com).

## Instalación

Para recrear el entorno de desarrollo, se pueden seguir los siguientes pasos:
#### Requisitos

 * RStudio
 * Python 3.9 o superior
 * Pip
    * Virtualenv (paquete de Python)
 * R 4.2 o superior
    * Renv (paquete de R)

#### Dependencias

En RStudio o en la terminal, al iniciar la consola interactiva de R dentro del directorio de la aplicación, se instalarán automáticamente las dependencias de R. Si esto no sucede, se puede utilizar el comando  `renv::restore()` dentro de la consola de R.

Como no se usó reticulate, las dependencias de Python deberán instalarse por separado. Se recomienda el uso de un entorno virtual de Python, pero se pueden instalar independientemente solo con pip. Para crear el entorno, se puede utilizar el siguiente comando:

```python3 -m venv env ```

Para ingresar al entorno, se puede utilizar el siguiente comando:

```source env/bin/activate ```

Para instalar las dependencias, se debe utilizar el siguiente comando:

``` pip3 install3 install -r requirements.txt ```

Para salir del entorno virtual, se puede utilizar el siguiente comando:
```deactivate```

La lista de paquetes utilizados y sus versiones en Python se pueden leer en el documento `requirements.txt` y en R en `renv.lock`.

Si alguno de los paquetes de R ha sido retirado de los repositorios de CRAN, se puede compilar desde su código fuente en GitHub utilizando los paquetes `remotes` o `devtools`.


## Descripción de los archivos del repositorio

Los descriptores se calcularon en los archivos `.py ` utilizando las moléculas que se encuentran en `SMILES.csv`. Los resultados de los descriptores se encuentran guardados en el directorio `DescriptorsCSV`. 

El código con el que se entrenó y validó los modelos se encuentra en los archivos `.R` identificados según su bioactividad. Es importante considerar que tomaron como entrada los archivos CSV que se encuentran en el directorio `DescriptorsCSV` (con los descriptores moleculares) y en el directorio `target` (con los vectores objetivos de cada bioactividad). La salida de datos de los modelos se guardo en los archivos en los archivos `.txt` que tienen como prefijo `output`.

Cada archivo `.R` de los modelos genera 3 objetos de R guardados en los archivos `.rds`. Los archivos que inician con `modelo` contienen el modelo final generado, los que empiezan con `predictores_filtrados` contienen las variables elegidas por ingeniería de características, y los que inician con `trained_recipe` contienen la receta del paquete `recipe` para el preprocesamiento de los datos con centrado y escalado.

Estos objetos y las salidas de datos usaron las iniciales de la bioactividad como nomenclatura de sufijos, siendo `ac` para anticancerígenos, `ad` para antidiabéticos, `ai` para antiinflamatorios, `am` para antimicrobianos y `ao` para antioxidantes.
