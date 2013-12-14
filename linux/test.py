#!/usr/bin/env 

from VectorMath import Vector

va = Vector((3.0, 4.0, 0.0))
vb = Vector((9.0, 2.1, 5.2))

print(va.components, " + ", vb.components)
vc = va + vb
print(vc.components)
print(va.components, " - ", vb.components)
vc = va - vb
print(vc.components)
print(va.components, " dot ", vb.components)
vc = va.dot(vb)
print(vc)
print(va.components, " x ", vb.components)
vc = va.cross(vb)
print(vc.components)
print(va.components, " ang ", vb.components)
vc = va.angle(vb)
print(vc)
print(va.components, " length ")
vc = va.length()
print(vc)
print(vb.components, " length ")
vc = vb.length()
print(vc)
print(va.components, " norm ")
vc = va.norm()
print(vc.components)
print(vb.components, " norm ")
vc = vb.norm()
print(vc.components)
