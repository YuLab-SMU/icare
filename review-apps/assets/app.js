async function initialiseApp() {
  const id = document.body.dataset.module;
  if (!id) return;
  const status = document.querySelector('#app-status');
  try {
    const response = await fetch('../assets/apps.json', {cache: 'no-store'});
    if (!response.ok) throw new Error('Configuration could not be loaded.');
    const config = await response.json();
    if (!config[id]) return;
    const url = new URL(config[id]);
    if (url.protocol !== 'https:' || !url.hostname.endsWith('.shinyapps.io') || url.username || url.password) {
      throw new Error('The application address is not a valid shinyapps.io HTTPS URL.');
    }
    const launch = document.querySelector('#launch');
    launch.href = url.href;
    launch.removeAttribute('aria-disabled');
    launch.removeAttribute('tabindex');
    status.textContent = 'Hosted on shinyapps.io · first launch may take a moment';
    document.querySelector('.dot').classList.add('ready');
    const placeholder = document.querySelector('#placeholder');
    placeholder.replaceChildren();
    const message = document.createElement('p');
    message.textContent = 'Open the application in a new tab, or load it below to keep the instructions nearby.';
    const embed = document.createElement('button');
    embed.className = 'button secondary';
    embed.textContent = 'Load application here';
    embed.addEventListener('click', () => {
      const frame = document.querySelector('#app-frame');
      frame.addEventListener('load', () => {
        status.textContent = 'Embedded application · use “Open application” if this view is unavailable';
      }, {once: true});
      frame.src = url.href;
      frame.hidden = false;
      placeholder.hidden = true;
      status.textContent = 'Application loading · use “Open application” if the embedded view is unavailable';
    }, {once: true});
    placeholder.append(message, embed);
  } catch (error) {
    status.textContent = 'Application configuration unavailable';
    document.querySelector('#placeholder p').textContent = error.message;
  }
}
initialiseApp();
