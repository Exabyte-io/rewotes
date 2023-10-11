import assert from "assert"
import React from 'react'
import { wrapInTestContext } from 'react-dnd-test-utils'
import { DragDropContext } from 'react-dnd'
import TestUtils from 'react-dom/test-utils'
import expect from 'expect'
import { App } from '/imports/ui/App.jsx'
import { DragComponent } from '/imports/ui/components/DragComponent'

describe("simple-todos-react", function () {
  it("package.json has correct name", async function () {
    const { name } = await import("../package.json")
    assert.strictEqual(name, "simple-todos-react")
  })

  if (Meteor.isClient) {
    it("client is not server", function () {
      assert.strictEqual(Meteor.isServer, false)
    })
  }

  if (Meteor.isServer) {
    it("server is not client", function () {
      assert.strictEqual(Meteor.isClient, false)
    })
  }
})

describe("drag-and-drop", function () {
  it("simple drag and drop", () => {
    // Render with the test context that uses the test backend
    const [AppContext, getBackend] = wrapInTestContext(App)
    const root = TestUtils.renderIntoDocument(<AppContext name="test" />)
  
    // Find the drag source ID and use it to simulate the dragging operation
    const dragComp = TestUtils.findRenderedDOMComponentWithTag(root, DragComponent)
    getBackend().simulateBeginDrag([dragComp.getHandlerId()])
  
  
    // See other TestBackend.simulate* methods for more!
  })
})
