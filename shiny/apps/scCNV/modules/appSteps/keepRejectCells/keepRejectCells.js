var cellToggleOverride = function(prefix, cell_id, key, override){
    Shiny.setInputValue(
        prefix + 'cellToggleOverride',
        {cell_id: cell_id, key: key, override: override},
        {priority: "event"}
    );
};
// var cellSetModalCN = function(prefix, cell_id, inputId){
//     Shiny.setInputValue(
//         prefix + 'cellSetModalCN',
//         {cell_id: cell_id, value: $("#" + inputId).val()},
//         {priority: "event"}
//     );
// };
