// enable a larger target area for track item settings icons
Shiny.addCustomMessageHandler('updateTrackDialogs', function(trackId) {
    $(".browserTrackItems, .browserTrackSettings").off("click").on("click", function(){
        const id = $(this).find("a").attr('id');
        Shiny.setInputValue(id, 1, {priority: "event"});
    });
});

// handle remove item clicks in trackItemsDialog
Shiny.addCustomMessageHandler('updateRemoveItemLinks', function(removeId) {
    $(".trackItemsTable a").off("click").on("click", function(){
        const tr = $(this).closest("tr");
        const key = $(tr).data("key");
        $(tr).remove();
        Shiny.setInputValue(removeId, key, {priority: "event"});
    });
});

// handle item option updates in trackItemsDialog
Shiny.addCustomMessageHandler('updateItemOptionChanges', function(config) {
    for(option in config.options){
        $(".trackItemsTable .changeItemOption").off("change").on("change", function(event){
            const tr = $(this).closest("tr");
            const key = $(tr).data("key");
            Shiny.setInputValue(
                config.id, 
                {
                    key: key,
                    option: option,
                    value: event.target.value
                }, 
                {priority: "event"}
            );
        });        
    }
});
