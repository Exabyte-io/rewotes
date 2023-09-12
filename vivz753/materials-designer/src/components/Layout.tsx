import React, { FC, useEffect, useState } from "react"
import Head from "next/head"
import { Explorer, SourceEditor, Toolbar, Visualizer } from "@components"

export const Layout: FC<React.PropsWithChildren> = ({ children }) => {
  const [input, setInput] = useState<string>("")

  useEffect(() => {
    console.log(input)
  }, [input])

  const handleEditor = async (input: string) => {
    setInput(input)
  }

  return (
    <>
      <Head>
        <title>Materials Designer</title>
        <meta name="description" content="prototype" />
        <meta name="viewport" content="width=device-width, initial-scale=1" />
        <link rel="icon" href="/favicon.ico" />
      </Head>
      <main>
        <div className="relative flex w-full flex-col overflow-auto bg-red-100">
          <Toolbar />
          <div className="flex min-h-screen w-full bg-blue-100 pt-20">
            <div className="flex w-full flex-row">
              <Explorer />
              <SourceEditor handleEditor={handleEditor} />
              <Visualizer input={input}/>
            </div>
            {/* {children} */}
          </div>
        </div>
      </main>
    </>
  )
}
