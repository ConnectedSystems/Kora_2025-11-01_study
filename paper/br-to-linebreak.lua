-- Convert HTML <br> / <br/> / <br /> raw inline elements to a Pandoc LineBreak
-- so they render correctly in docx (and other non-HTML) output formats.
function RawInline(el)
  if el.format == "html" and el.text:match("^<br") then
    return pandoc.LineBreak()
  end
end
