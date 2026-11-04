document.querySelectorAll('[data-copy]').forEach(button => {
  button.addEventListener('click', async () => {
    const code = document.getElementById(button.dataset.copy);
    const status = document.getElementById('copy-status');
    try {
      await navigator.clipboard.writeText(code.textContent);
      button.textContent = 'Copied';
      status.textContent = 'Code copied to clipboard';
    } catch {
      const selection = window.getSelection();
      const range = document.createRange();
      range.selectNodeContents(code);
      selection.removeAllRanges();
      selection.addRange(range);
      button.textContent = 'Copy manually';
      status.textContent = 'Code selected. Use your system copy shortcut.';
    }
    setTimeout(() => { button.textContent = 'Copy code'; }, 2200);
  });
});
if ('IntersectionObserver' in window) {
  const links = Array.from(document.querySelectorAll('.sidebar nav a'));
  const observer = new IntersectionObserver(entries => {
    const visible = entries.filter(entry => entry.isIntersecting);
    if (!visible.length) return;
    const id = visible[0].target.id;
    links.forEach(link => {
      if (link.hash === '#' + id) link.setAttribute('aria-current', 'true');
      else link.removeAttribute('aria-current');
    });
  }, {rootMargin: '-10% 0px -70% 0px'});
  document.querySelectorAll('.chapter').forEach(section => observer.observe(section));
}
const expandedForPrint = [];
window.addEventListener('beforeprint', () => {
  document.querySelectorAll('details:not([open])').forEach(item => {
    expandedForPrint.push(item); item.open = true;
  });
});
window.addEventListener('afterprint', () => {
  expandedForPrint.splice(0).forEach(item => { item.open = false; });
});
