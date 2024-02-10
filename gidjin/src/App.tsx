import { useState } from 'react'
import './App.css'

function App() {
  const [todo, setTodo] = useState([] as string[])
  const [input, setInput] = useState("")

  const addTodo = () => {
    setTodo([...todo, input])
    setInput("")
  }

  const onChange = (e: {target: {value: string}}) => {
    setInput(e.target.value)
  }

  const onDelete = (idxToRemove: number) => {
    setTodo(todo.filter((_item, idx) => idx !== idxToRemove))
  }

  return (
    <>
      <h1 className="text-3xl font-bold underline">
        Simple TODO app for demonstration purposes
      </h1>
      <ul className="pb-5">
        {
          todo.map((item, idx) => {
            return (
              <li data-testid={`item-${idx}`} key={idx}>
                {idx}. {item}
                <button data-testid={`delete-button-${idx}`} onClick={() => onDelete(idx)}>Delete</button>
              </li>)
          })
        }
      </ul>

      <input data-testid="add-input" type="text" onChange={onChange}value={input}/>
      <button data-testid="add-button" onClick={addTodo}>Add</button>

    </>
  )
}

export default App
