const Notes = () => {
	return (
		<div className='flex flex-col flex-nowrap gap-2 bg-slate-600 p-2 border-l-4 border-l-blue-500'>
			<p>
				<span className='font-bold italic'>Notes</span> - as you use the
				flowchart make sure you consider the following:
			</p>
			<ul className='list-disc list-inside'>
				<li>
					Nodes can only connect to one node at a time. Except for a conditional
					node.
				</li>
				<li>
					Conditional Nodes can go out to two outputs, depending on the
					truthState of the condition.
				</li>
				<li>
					The workflow will not execute unless the startNode and the endNode are
					connected through any arithmetic nodes.
				</li>
				<li>
					The test execution result will update when the workflow gets executed.
				</li>
			</ul>
		</div>
	);
};

export default Notes;
