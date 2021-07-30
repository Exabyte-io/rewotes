import React from "react"
import DragContainer from './DragContainer/DragContainer'
import DragComponent from './DragComponent/DragComponent'
import ElementContainer from './ElementContainer/ElementContainer'
import ComponentVisualizer from './ComponentVisualizer/ComponentVisualizer'
import "./Windows.scss"

/**
 * Component holding Three windows, 
 * The drag and drop options
 * The editable code
 * The visualizer during execution
 *
 * @param {*} {data, setData, runData, iterationIndex} data for code, function to edit data, data during execution, current index being executed
 */
const Windows = ({data, setData, runData, iterationIndex}) => {
  return (<div className="windows">
    <ElementContainer className={"drag-container"}>
      <DragComponent item={{type: "Declare", variable: "", value: ""}}>
        <div>{"Declare"}</div>
      </DragComponent>
      <DragComponent item={{type: "Add", variable: "", value: ""}}>
        <div>{"Add"}</div>
      </DragComponent>
      <DragComponent item={{type: "Subtract", variable: "", value: ""}}>
        <div>{"Subtract"}</div>
      </DragComponent>
      <DragComponent item={{type: "Multiply", variable: "", value: ""}}>
        <div>{"Multiply"}</div>
      </DragComponent>
      <DragComponent item={{type: "Divide", variable: "", value: ""}}>
        <div>{"Divide"}</div>
      </DragComponent>
      <DragComponent item={{type: "If", value1: "", value2: ""}}>
        <div>{"If"}</div>
      </DragComponent>
    </ElementContainer>
    <DragContainer data={data} setData={setData} iterationIndex={iterationIndex}/>
    <ComponentVisualizer state={runData}/>
  </div>)
}

export default Windows