import React from "react"
import "./SpeedButton.scss"

const SpeedButton = ({setSpeed, speed, children}) => {
  return (
    <div className="speed-button" onClick={() => setSpeed(speed)}>
      {children}
    </div>
  )
}

export default SpeedButton