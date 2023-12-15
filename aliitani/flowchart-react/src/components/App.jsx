import JsonViewer from './jsonViewer/JsonViewer';
import Flowchart from './flowchart/Flowchart';

import TestExecution from './TestExecution';
import WorkflowActions from './flowchart/WorkflowActions';

function App() {
	return (
		<div className='flex flex-nowrap flex-col justify-center gap-2 gap-x-5 items-center m-0 p-0 w-full'>
			<div className='flex flex-nowrap bg-slate-600 w-full px-2'>
				<p>
					Flowchart Diagram{' '}
					<small className='italic font-bold'>
						By Ali Itani (aliitani0@gmail.com)
					</small>
				</p>
			</div>
			<div className='flex flex-row flex-nowrap gap-x-2 gap-y-0'>
				<div className='flex flex-col flex-nowrap'>
					<JsonViewer />
					<TestExecution />
					<WorkflowActions />
				</div>
				<Flowchart />
			</div>
		</div>
	);
}

export default App;
