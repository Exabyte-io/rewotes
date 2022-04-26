import React from "react";
import Process from "./Process";

type Props = {
    id: string,
    data: string,
}

const Multiply: React.FC<Props> = (props: Props) => {
    const {id, data} = props


    return <>
        <Process isOperator={true} id={id} data={data} operator={'Multiply = '}/>
    </>
}

export default Multiply