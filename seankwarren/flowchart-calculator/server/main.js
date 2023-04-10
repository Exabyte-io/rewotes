import { Meteor } from 'meteor/meteor';
import { FlowsCollection } from '/imports/api/flows';

Meteor.startup(async () => {
    // We publish the entire LFlowsinks collection to all clients.
    // In order to be fetched in real-time to the clients
    Meteor.publish("flows", function () {
        return FlowsCollection.find();
    });
});

Meteor.methods({

    saveFlow(data) {
      FlowsCollection.insert(data);
    },

    fetchFlows() {
        return FlowsCollection.find().fetch();
    },
});
