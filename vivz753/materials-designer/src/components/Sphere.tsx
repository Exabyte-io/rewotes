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
  // const [active, setActive] = useState(false)

  // Subscribe this component to the render-loop, rotate the mesh every frame
  // useFrame(() => (mesh.current.rotation.x += 0.01))

  // const vertices = new Float32Array([
  //   -1.0,
  //   -1.0,
  //   1.0, // v0
  //   1.0,
  //   -1.0,
  //   1.0, // v1
  //   1.0,
  //   1.0,
  //   1.0, // v2

  //   1.0,
  //   1.0,
  //   1.0, // v3
  //   -1.0,
  //   1.0,
  //   1.0, // v4
  //   -1.0,
  //   -1.0,
  //   1.0, // v5
  // ])

  // Return view, these are regular three.js elements expressed in JSX
  return (
    <>
      <mesh
        {...props}
        ref={mesh}
        // scale={active ? 1.5 : 1}
        // onClick={() => setActive(!active)}
        onPointerOver={() => setHover(true)}
        onPointerOut={() => setHover(false)}
      >
        {/* <bufferGeometry>
        <bufferAttribute attach="attributes-position" count={vertices.length / 3}  array={vertices} itemSize={3} />
      </bufferGeometry> */}
        {/* <boxGeometry args={[1, 1, 1]} /> */}
        <sphereGeometry args={[0.1, 32, 16]} />
        <meshStandardMaterial color={hovered ? "white" : props.color} />
      </mesh>
    </>
  )
}
