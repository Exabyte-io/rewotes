import {NodeTypes} from 'react-flow-renderer';
import Terminal from "./Terminal";
import IO from "./IO";
import Process from "./Process";
import Decision from './Decision'
import Plus from "./Plus";
import Minus from "./Minus";
import Divide from "./Divide";
import Multiply from "./Multiply";

export type NodeNameTypes = 'terminal' | 'io' | 'process' | 'decision' | 'plus' | 'minus' | 'divide' | 'multiply'

const nodeTypes: NodeTypes = {
    terminal: Terminal,
    io: IO,
    process: Process,
    decision: Decision,
    plus: Plus,
    minus: Minus,
    divide: Divide,
    multiply: Multiply
};

export default nodeTypes