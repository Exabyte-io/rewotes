import { Toolbar } from "@components"
import Head from "next/head"
import React, { FC } from "react"

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
          {children}
        </div>
      </main>
    </>
  )
}
