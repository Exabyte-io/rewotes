import React from "react"
import classNames from "classnames"
import { v4 as uuid } from "uuid"
import "./DragContainer.scss"
import { DragContainerSlot } from "../DragContainerSlot"
import { isMathOp } from "../../utils"
import { MathComponent, DeclareComponent, IfComponent } from "../OperationSelections"

/**
 * The DragContainer renders operations from the json data.
 * The DragContainer also adds components to the stored json as appropriate when an element is dropped
 * 
 * @param {*} { className, data, setData } Additional classnames, the json data to be rendered, function to manipulate base json
 */
const DragContainer = ({ className, data, setData }) => {
  //Array containing unique keys of items stored in blocks
  const blockOrder = data.blockOrder
  //Object with key value pairs, key as the unique key of item, and item as item object
  const blocks = data.blocks

  /**
   * Adds item to data
   * Adds key of item to blockOrder at location index
   * Adds item to blocks with the key as it's key and value as item
   *
   * @param {*} key Key of the item
   * @param {*} index Location of item within the data
   * @param {*} item Item to be added
   */
  const addToBlocks = (key, index, item) => {
    setData( prevState => {
      return {
        blockOrder:  [
          ...prevState.blockOrder.slice(0, index),
          key,
          ...prevState.blockOrder.slice(index)
        ],
        blocks: {
          ...prevState.blocks,
          [key]: item
        }
      }
    })
  }

  /**
   * Edits the item with key as key, changes property name to value
   *
   * @param {*} key Key of the item
   * @param {*} name Property to be changed
   * @param {*} value Value property is to be changed to
   */
  const editBlock = (key, name, value) => {
    setData(prevState => {
      return {
        blockOrder: prevState.blockOrder, 
        blocks: {...prevState.blocks, [key]: {...prevState.blocks[key], [name]: value}}}
    })
  }

  const createSlot = (index, key) => {
    return <DragContainerSlot index={index} key={key} addToBlocks={addToBlocks}/>
  }

  /**
   * Renders all OperationSelections and DragContainerSlots
   *
   * @param {*} componentList Items 
   * @returns
   */
  const createComponentsfromList = (componentList, itemsObject) => {
    const retList = []
    
    componentList.forEach( (key, index) => {
      const curObj = itemsObject[key]
      if (curObj.type === "Slot") {
        retList.push(createSlot(index, curObj.key))
      }
      else if (isMathOp(curObj.type)) {
        retList.push(
          <MathComponent
            blockItem={curObj}
            editBlock={(name, value) => editBlock(key, name, value)}
            key={key}
          />
        )
      }
      else if (curObj.type === "Declare") {
        retList.push(
          <DeclareComponent
            blockItem={curObj}
            editBlock={(name, value) => editBlock(key, name, value)}
            key={key}
          />
        )
      }
      else if (curObj.type === "If") {
        retList.push(
          <IfComponent
            blockItem={curObj}
            editBlock={(name, value) => editBlock(key, name, value)}
            key={key}
          />
        )
      }
    })
    return retList
  }
  return (
    <div className={classNames("drag-container")}>
      {createComponentsfromList(blockOrder, blocks)}
    </div>
  )
}

export default DragContainer