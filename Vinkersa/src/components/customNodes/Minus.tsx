import React from "react";
import Process from "./Process";

type Props = {
    id: string,
    data: string,
}

const Minus: React.FC<Props> = (props: Props) => {
    const {id, data} = props


    return <>
        <Process isOperator={true} id={id} data={data} operator={'Minus = '}/>
    </>
}

export default Minus