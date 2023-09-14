"use client"
import { Line, Sphere } from "@components"
import { TrackballControls } from "@react-three/drei"
import { Canvas } from "@react-three/fiber"
import { FC } from "react"
import { Vector3 } from "three"

const colors = ["red", "yellow", "green", "purple", "blue"]

export const Visualizer: FC<{ latticeVectors: Vector3[]; pointVectors: Vector3[] }> = ({
  latticeVectors,
  pointVectors,
}) => {
  return (
    <div className="relative flex h-full w-full flex-col justify-center rounded-sm border border-dark1">
      <Canvas className="h-full w-full bg-dark2">
        <TrackballControls />
        <ambientLight />
        <pointLight position={[10, 10, 10]} />
        {pointVectors.map((vector, i) => {
          const index = i % colors.length
          return <Sphere key={i} position={vector} color={colors[index]} />
        })}
        <Line vectors={latticeVectors} />
      </Canvas>
      <p className="absolute top-0 mt-10 bg-dark2 px-10 text-center font-mozart text-xl uppercase tracking-widest text-light">
        Visualizer
      </p>
    </div>
  )
}
