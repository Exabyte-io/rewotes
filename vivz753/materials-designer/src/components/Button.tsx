"use client"

import { FC, ReactNode, MouseEventHandler } from "react"
import Link from "next/link"
// import clsx from "clsx"
import { ConditionalWrapper } from "@components"

export const Button: FC<{
  children: ReactNode
  href?: string
  size?: "xs" | "sm" | "lg"
  wide?: boolean
  className?: string
  variant?: "solid1"
  // onClick?: MouseEventHandler<HTMLButtonElement>
  disabled?: boolean
}> = ({ wide, children, href, size = "xs", className, variant = "solid1", disabled }) => {
  const sizeStyle =
    size === "xs" ? "p-2 lg:px-5" : size === "sm" ? "h-[49px] text-[20px] px-7" : "h-[66px] text-[24px] px-[36px]"
  const baseStyle = "select-none whitespace-nowrap rounded-md smooth-transition"
  const disabledStyle =
    "relative cursor-not-allowed after:bg-white after:w-full after:absolute after:h-full after:opacity-50 after:rounded-md"
  return (
    // <ConditionalWrapper condition={!!href} wrapper={(children) => <Link href={href ?? ""}>{children}</Link>}>
      <button
        disabled={disabled}
        // onClick={onClick}
        // className={clsx(
        //   "flex shrink-0 items-center justify-center",
        //   wide ? "w-full" : "max-w-max",
        //   baseStyle,
        //   sizeStyle,
        //   disabled && disabledStyle,
        //   variant === "solid1" && "bg-black text-white ring ring-orange-500",
        //   className
        // )}
      >
        {children}
      </button>
    // </ConditionalWrapper>
  )
}
