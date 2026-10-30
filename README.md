# rare-disease-simulator

Genera pacientes sintéticos con enfermedades poco frecuentes. A partir de las anotaciones de HPO arma el perfil de cada enfermedad y simula casos con sexo, edad y síntomas presentes y ausentes. Lo usamos en las primeras fases del reentrenamiento del modelo de [diagnostic.ar](https://github.com/MateAlv/diagnosticar-plataforma) y para evaluar el [mecanismo de repreguntas](https://github.com/MarBL00126/diagnosticar-question-agent).

Los datos de entrada (HPO, Orphanet y los perfiles de genes) son las copias fijadas en [diagnostic.ar-training](https://github.com/MateAlv/diagnostic.ar-training). El diseño y cada decisión están en [docs/README.md](docs/README.md), y los cambios de cada versión en [CHANGELOG.md](CHANGELOG.md).

Es parte del Trabajo Profesional de Ingeniería en Informática (FIUBA) de Mateo Alvarez, Mariano Brizuela López, Gonzalo Bustamante y Lucía Lourengo Caridade.

## Probarlo

```bash
pip install -e ".[dev]"
pytest
```

Para generar casos hace falta `diagnostic.ar-training` clonado al lado de este repo, porque de ahí salen los datos y las calibraciones. Los comandos y sus opciones están en [docs/README.md](docs/README.md).

## Datos

Usa la ontología HPO y sus anotaciones, y datos de Orphanet (CC BY 4.0). Todos los casos son sintéticos: no hay datos de pacientes reales.

## Licencia

Apache-2.0.
