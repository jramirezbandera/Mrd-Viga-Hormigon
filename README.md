# MRd Fibers EC2

Script en Python para calcular el momento resistente de cálculo (MRd) de secciones rectangulares de hormigón armado mediante el método de fibras, según el Eurocódigo 2 (EC-2).

## Características
- Cálculo por integración de fibras de hormigón y capas de acero.
- Compatible con hormigones de fck ≤ 50 MPa.
- Implementa diagrama parabólico-rectangular para hormigón y diagrama elasto-plástico para acero.
- Permite calcular MRd positivo y negativo.

## Requisitos
- Python 3.8 o superior.
- Librerías estándar: `math`, `dataclasses`, `typing`.

## Uso
Ejecutar directamente el script:

```bash
python Mrd.py
```

Ejemplo de salida:
```
MRd+ = 307.85 kN·m (c=42.3 mm, z≈545.7 mm)
MRd- = 308.02 kN·m (c=42.4 mm, z≈545.6 mm)
```

## Licencia
Este proyecto se distribuye bajo la licencia MIT. Consulta el archivo [LICENSE](LICENSE) para más información.
