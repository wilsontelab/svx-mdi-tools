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
// initialize plot interactions

Shiny.addCustomMessageHandler("cellPlotsWrapperInit", function(opt){

    // parse interactive targets 
    const wrapperId = "#" + opt.prefix + "cellPlotsWrapper";
    const mdiVertical = wrapperId + " .cellStackVertical";

    // convert document to widget coordinates
    const relCoord = function(event, element){return {
        x: event.pageX - $(element).offset().left,
        y: event.pageY - $(element).offset().top,
    }}

    // handle all requested interactions, listed here if rough order of occurrence
    $(wrapperId).off("mouseenter").on("mouseenter", function() {
        $(mdiVertical).show();
    });
    $(wrapperId).off("mousemove").on("mousemove", function(event) {
        const coord = relCoord(event, this);
        $(mdiVertical).css({left: coord.x - 1});
    });
    $(wrapperId).off("mouseleave").on("mouseleave", function() {
        $(mdiVertical).hide();
    });
});
Shiny.addCustomMessageHandler("cellPlotsWrapperUpdate", function(opt){
    const wrapperId = "#" + opt.prefix + "cellPlotsWrapper";
    const mdiVertical = wrapperId + " .cellStackVertical";
    $(mdiVertical).css({height: (1.35 * 96 + 1) * opt.cellsPerPage});
});
