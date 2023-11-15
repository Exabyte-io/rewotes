import { MeshProps } from "@react-three/fiber"
import { useRef, useState } from "react"
import { Mesh } from "three"

type SphereProps = {
  color: string
}

export function Sphere(props: MeshProps & SphereProps) {
  // This reference will give us direct access to the mesh
  const mesh = useRef<Mesh>(null!)
  // Set up state for the hovered and active state
  const [hovered, setHover] = useState(false)

  return (
    <>
      <mesh {...props} ref={mesh} onPointerOver={() => setHover(true)} onPointerOut={() => setHover(false)}>
        <sphereGeometry args={[0.05, 32, 16]} />
        <meshStandardMaterial color={hovered ? "black" : props.color} />
      </mesh>
    </>
  )
}
