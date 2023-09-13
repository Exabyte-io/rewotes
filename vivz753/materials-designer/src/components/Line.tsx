import React, { FC, useLayoutEffect } from "react"
import { BufferGeometry, Vector3 } from "three"

export const Line: FC<{ vectors: Vector3[]; color?: string }> = ({ vectors, color }) => {
  const points: Vector3[] = vectors || [new Vector3(-1, 0, 0), new Vector3(0, 1, 0), new Vector3(1, 0, 0)]

  const ref = React.useRef<BufferGeometry>(null!)
  useLayoutEffect(() => {
    if (ref.current) {
      ref.current.setFromPoints(points)
    }
  }, [points])

  return (
    <line>
      <bufferGeometry attach="geometry" ref={ref} />
      <lineBasicMaterial color={color || "#9c88ff"} linewidth={10} linecap={"round"} linejoin={"round"} />
    </line>
  )
}
