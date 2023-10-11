import React, { useState } from "react"
import {v4 as uuid} from "uuid"
import "./App.scss"
import { DndProvider } from "react-dnd"
import { HTML5Backend } from "react-dnd-html5-backend"
import { RunButtons } from "./components/Header"
import Windows from "./components/Windows"
import { isNil } from "lodash"
import { checkIfNumber } from "./utils"
import { BrowserRouter as Router, Route } from 'react-router-dom'


const App = () => {
  //Initial Data for code
  const initialKey = uuid()
  const [data, setData] = useState({
    blockOrder: [initialKey],
    blocks: { [initialKey]: { "type": "Slot" }} })

  //Speed of code execution
  const [speed, setSpeed] = useState(1000)

  //Current location of code execution
  const [iterationIndex, setIterationIndex] = useState(0)

  //Json representation for the running code
  const [runData, setRunData] = useState({vars: {}})

  //Check if we have exited an if statement
  const [prevSlotted, setPrevSlotted] = useState(false)
  //If if statement and all components inside should be skipped
  const [ifsInvalid, setIfsInvalid] = useState(0)

  const timer = ms => new Promise(res => {
    setTimeout(res, ms)
  })

  /**
   * Resets code execution
   *
   */
  const reset = () => {
    setRunData({vars: {}})
    setIterationIndex(0)
  }


  /**
   * Executes the code according to data, our source of truth for 
   *
   */
  async function runCode () {
    for (const key of data.blockOrder) {
      const newData = parseLine(runData, data.blocks[key])
      if (newData === "Error") {
        setRunData("Error")
        setIterationIndex(0)
        break
      }
      else {
        setRunData({...newData})
        setIterationIndex(prev => prev + 1)
        await timer(speed)
      }
      
    }
    setIterationIndex(0)
  }

  /**
   * Traverses up the stack frame and sets the respective variable to the given value
   *
   * @param {*} curData Current value of the running code
   * @param {*} varName Variable name
   * @param {*} varValue value of the variable
   */
  const setVariableValue = (curData, varName, varValue) => {
    if (!isNil(curData)) {
      if (!isNil(curData.vars[varName])) {
        curData.vars[varName] = varValue
      }
      else {
        setVariableValue(curData.parentFrame, varName, String(varValue))
      }
    }
  }


  /**
   * Gets the value of the variable by traversing up the tree.
   * Takes literal value of variable name if the variable is not found in an earlier stackframe
   *
   * @param {*} curData Current value of the running code
   * @param {*} varName Variable name
   * @returns
   */
  const getVariableValue = (curData, varName) => {
    if (isNil(curData)) {
      return varName
    }
    
    if (!isNil(curData.vars[varName])) {
      return curData.vars[varName]
    }
    else {
      return getVariableValue(curData.parentFrame, varName)
    }
  }
    
  /**
   * Parses a single line of code by inspecting the item, returns the updated running state
   *
   * @param {*} curData Current state of the running application
   * @param {*} curItem Current line/item to be executed
   * @returns Updated State
   */
  const parseLine = (curData, curItem) => {
    if (curItem.type === "Slot") {
      if (prevSlotted) {
        setIfsInvalid(prev => Math.max(prev - 1, 0))
        return curData.parentFrame
      }
      //Check if we have exited an if statement via a double slot
      setPrevSlotted(true)
    }
    else {
      setPrevSlotted(false)
    }

    //if ifsInvalid is above 0, we are still in a for loop we must skip, though we do not change the state
    if (ifsInvalid > 0) {
      return curData
    }
    else if (curItem.type.match(/^(Add|Subtract|Multiply|Divide)$/)) {
      let varValue = getVariableValue(curData, curItem.variable)
      let value = getVariableValue(curData, curItem.value)

      if (!checkIfNumber(value) || !checkIfNumber(varValue)) {
        return "Error"
      }
      else {
        varValue = parseFloat(varValue)
        value = parseFloat(value)

        if (curItem.type === "Add") {
          setVariableValue(curData, curItem.variable, varValue + value)
        }
        else if (curItem.type === "Subtract") {
          setVariableValue(curData, curItem.variable, varValue - value)
        }
        else if (curItem.type === "Multiply") {
          setVariableValue(curData, curItem.variable, varValue * value)
        }
        else if (curItem.type === "Divide") {
          if (value === "0") {
            return "Error"
          }
          setVariableValue(curData, curItem.variable, varValue / value)
        }
      }
      
    }
    else if (curItem.type === "Declare") {
      if (!isNil(curData.vars[curItem.variable])) {
        return "Error"
      }
      curData.vars[curItem.variable] = getVariableValue(curData, curItem.value)
    }
    else if (curItem.type === "If") {
      if (getVariableValue(curData, curItem.value1) === getVariableValue(curData, curItem.value2)) {
        curData = {
          vars: {},
          parentFrame: {...curData}
        }
      }
      else {
        setIfsInvalid(prev => prev + 1)
      }
    }
    return curData
  }

  return (
    <DndProvider backend={HTML5Backend}>
      <div className={"App"}>
        <RunButtons setSpeed={setSpeed} runCode={runCode} reset={reset}/>
        <Windows data={data} setData={setData} iterationIndex={iterationIndex} runData={runData}/>
      </div>
    </DndProvider> 
  )
}
  
export default App