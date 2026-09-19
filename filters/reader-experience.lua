-- Public reading experience; leaves calculations, outputs and figure alt text intact.
-- Quarto runs this filter before its built-in HTML rendering filters.
local function fold_code(el)
  if el.attributes and el.attributes['code-fold'] == 'show' then
    el.attributes['code-fold'] = 'true'
    return el
  end
end

function Pandoc(doc)
  if not FORMAT:match('html') then return doc end
  doc = doc:walk({Div = fold_code, CodeBlock = fold_code})

  local function panel(key, id, class)
    local value = doc.meta[key]
    if not value then return nil end
    local kind = pandoc.utils.type(value)
    local blocks
    if kind == 'Blocks' then
      blocks = value
    elseif kind == 'Inlines' then
      blocks = {pandoc.Para(value)}
    else
      blocks = pandoc.read(pandoc.utils.stringify(value), 'markdown').blocks
    end
    return pandoc.Div(blocks, pandoc.Attr(id, {class}))
  end

  local intro = panel('reader-intro', 'article-takeaways', 'reader-takeaways')
  if intro then doc.blocks:insert(1, intro) end

  local contact = panel('reader-contact', 'article-contact', 'reader-contact')
  if contact then
    local position = #doc.blocks + 1
    for i, block in ipairs(doc.blocks) do
      if block.t == 'Header' and block.identifier == 'references' then
        position = i
        break
      end
    end
    doc.blocks:insert(position, contact)
  end
  return doc
end
