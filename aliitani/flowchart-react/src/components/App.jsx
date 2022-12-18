import JsonViewer from './jsonViewer/JsonViewer';
import Flowchart from './flowchart/Flowchart';

import TestExecution from './TestExecution';

function App() {
	return (
		<div className='flex flex-nowrap flex-col justify-center gap-y-3.5 gap-x-5 items-center m-0 p-0 w-full'>
			<div className='flex flex-nowrap bg-slate-600 w-full px-2'>
				<p>
					Flowchart Diagram{' '}
					<small className='italic font-bold'>
						By Ali Itani (aliitani0@gmail.com)
					</small>
				</p>
			</div>
			<div className='flex flex-row flex-nowrap gap-x-2 my-5'>
				<div>
					<JsonViewer />
					<TestExecution />
				</div>
				<Flowchart />
			</div>
		</div>
	);
}

export default App;
