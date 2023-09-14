"use client"
import { Line, Sphere } from "@components"
import { TrackballControls } from "@react-three/drei"
import { Canvas } from "@react-three/fiber"
import { FC } from "react"
import { Vector3 } from "three"

const colors = ["red", "yellow", "green", "purple"]

export const Visualizer: FC<{ vectors: Vector3[] }> = ({ vectors }) => {
  return (
    <div className="flex h-full w-full flex-col justify-center">
      <p className="text-center text-xl">Visualizer</p>
      <Canvas className="h-full w-full bg-pink-200">
        <TrackballControls />
        <ambientLight />
        <pointLight position={[10, 10, 10]} />
        {vectors.map((vector, i) => (
          <Sphere key={i} position={vector} color={colors[i - 1 % (colors.length - 1)]} />
        ))}
        <Line vectors={vectors} />
      </Canvas>
    </div>
  )
}
