using System.Web.Mvc;

namespace PhdThesisIngaSamonenko.Controllers
{
    public class GenzController : Controller
    {
        public ActionResult Index()
        {
            return File("Genz.R", "text/plain");
        }

    }
}
