---
title: How To
layout: template
filename: test
--- 

# How to Create a Multi-page Website using Github Pages

### To Create the Template
1. Create the page normally following [pages.github.com](https://pages.github.com), but use `CONTENT` as the content of the page
2. Choose a theme and publish the page
3. Fetch and checkout the gh-pages branch on your local repository
4. Create a directory `_layouts` in the repository
5. Rename `index.html` to `template.html` and move it into the `_layouts` directory
6. Open `template.html` and replace the `<p>CONTENT</p>` placeholder with {% raw %}`{{ content }}`{% endraw %} (this is [Jekyll](https://jekyllrb.com) syntax to grab the content from the MarkDown pages you will create)
7. Identify the navigation/button section of HTML
8. Copy one navigation/button item (probably a `<a href="">` or similar tag)
9. Insert this code at the top of the navigation/button item section:

```
{% raw %}
{% for page in site.pages %}
    <a href={{ page.filename }}>{{ page.title }}</a>
{% endfor %}
{% endraw %}
```

To match your theme, paste the copied navigation/button item in place of `<a href={{ page.filename }}>{{ page.title }}</a>`, but use `{{ page.filename }}` for the href and `{{ page.title }}` for the content (as shown in the example in step 9)

### To Create Your First Page
1. Make a new file called `index.md` in your repository
2. Copy the content of your `readme.md` or write a new home page in MarkDown into this file
3. At the top of this file, add the following:

```
---
title: PAGE TITLE HERE
layout: template
filename: NAME OF THIS .md FILE HERE
--- 
```

Commit your changes and push them to the gh-pages branch

Now, when you go to `YOURGITHUBNAME.github.io/YOURPROJECTNAME`, you should see the contents of your index.md formatted with the theme that you chose.

### To Create Additional Pages
1. Make a new file called `PAGENAME.md` in your repository (where PAGENAME is the name of your new page)
2. Write the content for this new page in MarkDown
3. At the top of this file, add the following:

```
---
title: PAGE TITLE HERE
layout: template
filename: NAME OF THIS .md FILE HERE
--- 
```

Commit your changes and push them to the gh-pages branch

Now, when you go to `YOURGITHUBNAME.github.io/YOURPROJECTNAME`, you should see a link to your new page. When you click this link, you should see your new page formatted with the theme that you chose.