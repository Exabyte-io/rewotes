import React, { FC } from "react"
import Head from "next/head"

export const Layout: FC<React.PropsWithChildren> = ({ children }) => {
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
              <SourceEditor />
              <Visualizer />
            </div>
            {/* {children} */}
          </div>
        </div>
      </main>
    </>
  )
}

const Toolbar = () => {
  return (
    <div className="absolute top-0 z-10 flex h-20 max-h-20 w-full flex-row items-center justify-start gap-5 bg-green-200 px-20">
      Toolbar
    </div>
  )
}

const Explorer = () => {
  return (
    <div className="flex h-full max-w-max flex-col whitespace-nowrap bg-gray-500 p-10">
      <p className="text-xl">Explorer</p>
    </div>
  )
}

const SourceEditor = () => {
  return (
    <div className="flex h-full max-w-max flex-col whitespace-nowrap bg-gray-800 p-10">
      <p className="text-xl">Source Editor</p>
    </div>
  )
}

const Visualizer = () => {
  return <div className="h-full w-full grow bg-yellow-200"></div>
}
